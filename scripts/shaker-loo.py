#!/usr/bin/env python
import pickle
import argparse
from collections import defaultdict
import shaker.rna_tools.rna_io as rio
import shaker.simushape as sim
import shaker.rna_tools.util as util
from eden import graph as eg
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
import numpy as np
from utilities import loadRecords
import sys

parser = argparse.ArgumentParser(description='Perform LOO Cross-Validation with SHAKER')
parser.add_argument('--input', '-i',required=True, help="Input fasta like file")
parser.add_argument('--performance','-p',required=True,help="Result of cross validation")
parser.add_argument('--reactivity','-r',required=True,help="Predicted reactivities")
args = parser.parse_args()

def RMSE(predictions, targets):
    predictions = np.array(predictions)
    targets = np.array(targets)
    return np.sqrt(((predictions - targets) ** 2).mean())

def AUC(structure,reactivity):
    structure = np.array(list(structure))
    unpaired_mask = (structure == ".")
    unpaired = np.zeros(structure.shape[0])
    unpaired[unpaired_mask] = 1
    idx = np.where(~np.isnan(reactivity))[0]
    return roc_auc_score(unpaired[idx],reactivity[idx])

def main(args):
    print "Load input data ..."
    records = loadRecords(args.input) 
    data = {}
    for name in records.keys():
        data[name] = [records[name]["reactivity"],records[name]["sequence"],records[name]["structure"]]
        reactivity = []
        for x in data[name][0]:
            if np.isnan(x):
                reactivity.append(None)
            else:
                reactivity.append(x)
        data[name][0] = reactivity
    print "Done ."
    print "Train SHAKER model ..."
    
    fperformance = open(args.performance,"w") if args.performance != "-" else sys.stdout
    fperformance.write("\t".join(["name","spearmanr","p-value","AUROC-observed-reactivity","AUROC-predicted-reactivity","RMSE"])+"\n")
    
    fout = open(args.reactivity,"w")
    
    for name in data.keys():
        print name;
        keys = set(data.keys())
        keys.remove(name)
        # data[name][0] reactivity
        # data[name][1] sequence
        # data[name][2] structure
        model  = sim.make_model(data,list(keys))
        graph = util.sequence_dotbracket_to_graph(data[name][1],data[name][2])
        embedding = eg.vertex_vectorize([graph])[0]
        reactivity_pred = model.predict(embedding).reshape(-1)
        fout.write(">"+name+"\n")
        fout.write(",".join(np.round(reactivity_pred,3).astype(str))+"\n")
        reactivity = np.array(data[name][0]).astype(float)
        structure = data[name][2]
        auc = AUC(structure,reactivity)
        auc_pred = AUC(structure,reactivity_pred)
        nan_mask = np.isnan(reactivity)
        reactivity = reactivity[~nan_mask]
        reactivity_pred = reactivity_pred[~nan_mask]
        corr,p = spearmanr(reactivity_pred,reactivity)
        rmse = RMSE(reactivity_pred,reactivity) 
        fperformance.write("\t".join([name,str(corr),str(p),str(auc),str(auc_pred),str(rmse)])+"\n")
    fperformance.close()
    fout.close()



if __name__ == "__main__":
    main(args)

