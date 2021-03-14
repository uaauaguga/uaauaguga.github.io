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
from utilities import loadRecords,annotatePairing,getKmer,writeRecords


def prepareDataset(dataDict):
    dataset = {}
    for seq_id in dataDict.keys():
        structure = dataDict[seq_id]["pairing"]
        reactivity = dataDict[seq_id]["reactivity"]
        dataset[seq_id] = getKmer(structure,reactivity,flankingL = 2)
    return dataset

def fitting(dataset,frequency=200):
    reactivities = []
    
    # Summary instances into a table
    for seq_id in dataset.keys():
        structure_contexts,reactivity = dataset[seq_id]
        L = len(reactivity)
        assert len(structure_contexts) == L, "Structure context and reactivities has different length"
        for i in range(L):
            reactivities.append((structure_contexts[i],reactivity[i],seq_id))
    reactivities = pd.DataFrame.from_records(reactivities)
    reactivities.columns = ["structure-context","reactivity","sequence-id"]
    # Only use reactivity larger than zero
    reactivities = reactivities[reactivities["reactivity"]>0]
    # Get instances number of each structure context
    structure_contexts = reactivities["structure-context"].unique()
    n_context = structure_contexts.shape[0]
    n_instances = reactivities.groupby("structure-context").apply(lambda x:x.shape[0])
    statistics = pd.DataFrame(index=structure_contexts,columns=["instances","assignment"])    
    statistics.loc[n_instances.index,"instances"] = n_instances.values

    print("5 mer should have 32 structure contexts")
    frequent_set = set(n_instances[n_instances>frequency].index)
    print("{} present in the input dataset".format(n_context))
    print("{} meet the frequency cutoff".format(len(frequent_set)))
    # Split structure context to frequent ones and not frequent ones    
    frequent_reactivities = reactivities[reactivities["structure-context"].isin(frequent_set)]
    not_frequent_reactivities = reactivities[~reactivities["structure-context"].isin(frequent_set)]

    ## Fit an individual generalized extreme distribution for each frequent 5 mer instance
    ## For rare instance, assign the instance the most similar fitted instance (fitted model with highest likelihood)
    # Structure context to model mapping
    modelDict = {}
    # logged likelihod of rare structure context to assign to a structure context
    likelihoodsDict = {}
    # Structure context not fit well to generalize extreme value distribution
    dubious_instances = []

    for structure_context in frequent_reactivities["structure-context"].unique():
        react = reactivities[reactivities["structure-context"]==structure_context]["reactivity"].values      
        shape,location,scale = genextreme.fit(react)
        if shape > 0:
            dubious_instances.append(structure_context)
            continue
        model = genextreme(shape,location,scale)
        modelDict[structure_context] = model
        likelihoodsDict[structure_context] = not_frequent_reactivities.groupby("structure-context").apply(lambda x:np.log(model.pdf(x["reactivity"].values)).sum())

    # Summarize the likelihod of rare structure context into a DataFrame
    likelihoods = pd.DataFrame(likelihoodsDict)
    # Map rare structure context to most similar frequent structure context
    assignment0 = dict(likelihoods.idxmax(axis=1))

    # Map frequent structure context to list of rare structure contexts
    assignment = defaultdict(list)
    for less,more in assignment0.items():
        statistics.loc[less,"assignment"] = more
        assignment[more].append(less)

    # Refit genextreme model for frequent structure context
    for more,less in assignment.items():
        for context in less:
            modelDict[context] = modelDict[more]
        contexts = set(less)
        contexts.add(more)
        react = reactivities[reactivities["structure-context"].isin(contexts)]["reactivity"].values
        shape,location,scale = genextreme.fit(react)
        model = genextreme(shape,location,scale)
        for structure_context in contexts:
            modelDict[structure_context] = model
    if len(dubious_instances) > 0:
        print("Fitting for the following instances generates positive shape value, which is dubious")
        print(",".join(dubious_instances))
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
            print("{} is assigned to {}".format(structure_context,max_instance))
            statistics.loc[structure_context,"assignment"] = max_instance
            modelDict[structure_context] = modelDict[max_instance]   

    #reactivities.to_csv("reactivity-table.txt",index=False,sep="\t")
    
    return modelDict,statistics
        

def main():
    parser = argparse.ArgumentParser(description='Simulate reactivity from generalized extreme value distribution of different 5mers')
    parser.add_argument('--dataset', '-d',required=True, help="Input fasta like file contains sequence,structure and reactivity")
    parser.add_argument('--frequency', '-f',type=int,default=200,help="Merge structure context with instance lower than this value")
    parser.add_argument('--statistics', '-s',required=True,help="Statistics for instance assignment")
    parser.add_argument('--model','-m',required=True,help="Model output")
    args = parser.parse_args()

    print("Load input data ...")
    dataDict = loadRecords(args.dataset,order="sequence,structure,reactivity",dtype="str,str,float")
    dataDict = annotatePairing(dataDict) 
    
    print("Convert dataset into instances of 5 mers ...")
    dataset = prepareDataset(dataDict)

    print("Model fitting ...")
    modelDict,statistics = fitting(dataset,args.frequency)
    statistics.to_csv(args.statistics,sep="\t")
    
    with open(args.model, 'wb') as f:
        pickle.dump(modelDict, f)
 
    
if __name__ == "__main__":
    main()


