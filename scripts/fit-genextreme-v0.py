#!/usr/bin/env python
import argparse
import numpy as np
from scipy.stats import genextreme
from collections import defaultdict
import pandas as pd
import pickle
import re
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
from utilities import loadRecords

def getKmer(X,y,flankingL = 2):
    isnan = np.isnan(y).astype(int)
    window = flankingL*2 + 1
    # any_flanking_nan: if this is at least nan in [-2,3], annotate the position as 1
    any_flanking_nan = np.convolve(isnan,np.ones(window,dtype=int),'valid')
    any_flanking_nan = np.concatenate([np.array([window]*flankingL),any_flanking_nan,np.array([window]*flankingL)])
    mask = any_flanking_nan==0
    features = []
    for idx in np.where(mask)[0]:
        x = X[idx-flankingL:idx+flankingL+1].reshape(-1)
        features.append(x)
    X_ = np.array(features)
    y_ = y[mask]
    return X_,y_

def annotatePairing(dataDict):
    for seq_id in dataDict.keys():
        annotation = dataDict[seq_id]["structure"].replace(".","U")
        annotation = re.sub("[^U]","P",annotation)
        dataDict[seq_id]["pairing"] = annotation
    return dataDict

def prepareDataset(records):
    dataset = {}
    for seq_id in records.keys():
        structure = np.array(list(records[seq_id]["structure"]))
        X = np.zeros(structure.shape[0])
        unpaired = (structure==".")
        X[~unpaired] = 1
        y = records[seq_id]["reactivity"]
        dataset[seq_id] = getKmer(X,y,flankingL = 2)
    return dataset


def digit2string(X):
    s = np.empty_like(X).astype(str)
    s[X==1],s[X==0] = "P","U"
    return np.apply_along_axis(lambda x:"".join(x),arr=s,axis=1)



def fitting(dataset,frequency=200):
    reactivities = []
    for seq_id in dataset.keys():
        X,y = dataset[seq_id]
        structure_contexts = digit2string(X)
        for i in range(y.shape[0]):
            reactivities.append((structure_contexts[i],y[i],seq_id))
    reactivities = pd.DataFrame.from_records(reactivities)
    reactivities.columns = ["structure-context","reactivity","sequence-id"]
    
    n_instances = reactivities.groupby("structure-context").apply(lambda x:x.shape[0])
    frequent_set = set(n_instances[n_instances>200].index)
    
    frequent_reactivities = reactivities[reactivities["structure-context"].isin(frequent_set)]
    not_frequent_reactivities = reactivities[~reactivities["structure-context"].isin(frequent_set)]

    ## Fit an individual generalized extreme distribution for each frequent 5 mer instance
    ## For rare instance, assign the instance the most similar fitted instance (fitted model with highest likelihood)
    modelDict = {}
    likelihoodsDict = {}

    dubious_instances = []

    for structure_context in frequent_reactivities["structure-context"].unique():
        react = reactivities[reactivities["structure-context"]==structure_context]["reactivity"]
        shape,location,scale = genextreme.fit(react)
        if shape > 0:
            dubious_instances.append(structure_context)
            continue
        model = genextreme(shape,location,scale)
        modelDict[structure_context] = model
        likelihoodsDict[structure_context] = not_frequent_reactivities.groupby("structure-context").apply(lambda x:np.log(model.pdf(x["reactivity"].values)).sum())
    likelihoods = pd.DataFrame(likelihoodsDict)
    assignment0 = dict(likelihoods.idxmax(axis=1))
    assignment = defaultdict(list)

    for less,more in assignment0.items():
        assignment[more].append(less)

    for more,less in assignment.items():
        contexts = set(less)
        contexts.add(more)
        react = reactivities[reactivities["structure-context"].isin(contexts)]["reactivity"].values
        shape,location,scale = genextreme.fit(react)
        model = genextreme(shape,location,scale)
        for structure_context in contexts:
            modelDict[structure_context] = model

    if len(dubious_instances) > 0:
        print("#"*50)
        print("Fitting for the following instances generates positive shape value, which is infeasible")
        print(",".join(dubious_instances))
        print("#"*50)
        max_ll = np.nan
        for structure_context in dubious_instances:
            for fitted_context,model in modelDict.items():
                x = reactivities[reactivities["structure-context"]==fitted_context]["reactivity"].values
                current_ll = np.log(model.pdf(x)).sum()
                if np.isnan(max_ll):
                    max_ll = np.log(model.pdf(x)).sum()
                    max_instance = fitted_context
                else:
                    if max_ll < current_ll:
                        max_ll = current_ll
                        max_instance = fitted_context
            print("This instance is assigned to {}".format(max_instance))
            modelDict[structure_context] = modelDict[max_instance]   
    return modelDict
        

def simulate(pairing,model_dict):
    reactivity = np.zeros(len(pairing))
    pairing = "UU" + pairing + "UU"
    for idx in range(len(reactivity)):
        reactivity[idx] = model_dict[pairing[idx:idx+5]].rvs()
    reactivity[reactivity<0] = 0
    smoothed = np.convolve(reactivity,np.array([1/6,2/3,1/6]), mode='valid')
    reactivity = np.concatenate([np.array([reactivity[0]]),smoothed,np.array([reactivity[-1]])])
    return reactivity
        
    


def main():
    parser = argparse.ArgumentParser(description='Simulate reactivity from generalized extreme value distribution of different 5mers')
    parser.add_argument('--input', '-i',required=True, help="Input fasta like file contains sequence,structure and reactivity")
    parser.add_argument('--performance','-p',required=True,help = "Performance table")
    parser.add_argument('--reactivity','-r',required=True, help = "Simulated reactivity")
    parser.add_argument('--model','-m',required=True,help="Model output")
    args = parser.parse_args()

    print("Load input data ...")
    records = loadRecords(args.input,data_type ="sequence,structure,reactivity")
    print("Done .")
    
    print("Convert dataset into instances of 5 mers ...")
    dataset = prepareDataset(records)
    print("Done .")
    
    print("Model fitting ...")
    model_dict = fitting(dataset)
    print("Done .")

    records = annotatePairing(records)
    performances = []
    
    freact = open(args.reactivity,"w")   

    for seq_id in records.keys():
        pairing = records[seq_id]["pairing"]
        predicted = simulate(pairing,model_dict)
        freact.write(">" + seq_id + "\n")
        freact.write(",".join(np.round(predicted,3).astype(str))+"\n")
        y = records[seq_id]["reactivity"]
        paired = np.array(list(pairing))
        paired[paired=="P"] = 1
        paired[paired=="U"] = 0
        paired = paired.astype(int)
        mask = np.isnan(y)
        y = y[~mask]
        y_pred = predicted[~mask]
        paired = paired[~mask]
        r = spearmanr(y,y_pred)[0]
        auroc_true = roc_auc_score(1-paired,y)
        auroc_pred = roc_auc_score(1-paired,y_pred)
        performances.append((seq_id,r,auroc_true,auroc_pred))
    performances = pd.DataFrame.from_records(performances)
    performances.columns = ["RNA","spearmanr","AUROC","AUROC-predicted"]
    performances = performances.set_index("RNA")
    performances.to_csv(args.performance,sep="\t")
    
    freact.close()
    print("Save fitted model ...")
    with open(args.model, 'wb') as f:
        pickle.dump(model_dict, f)
    print("Done .")

            
    

if __name__ == "__main__":
    main()


