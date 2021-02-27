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

import sys

parser = argparse.ArgumentParser(description='Train SHAKER model')
parser.add_argument('--input', '-i',required=True, help="Input fasta like file")
parser.add_argument('--output','-o',required=True,help="Result of cross validation")
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


def load_records(path,data_type ="sequence,structure,reactivity"):
    data_types = data_type.split(",")
    n_data_types = len(data_types)
    records = defaultdict(dict)
    with open(path) as f:
        for line in f:
            try:
                assert line.startswith(">")
            except:
                print "The input data is not consistent to specified data_type parameter"
            seq_id = line.strip()[1:].strip()
            for data_type in data_types:
                data = next(f).strip()
                if data_type == "reactivity":
                    data = np.array(data.split(",")).astype(float)
                records[seq_id][data_type] = data
    return records


def main(args):
    print "Load input data ..."
    records = load_records(args.input) 
    data = {}
    for name in records.keys():
        print records[name].keys()
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
    
    fout = open(args.output,"w") if args.output != "-" else sys.stdout
    fout.write("\t".join(["name","spearmanr","p-value","AUROC-observed-reactivity","AUROC-predicted-reactivity","RMSE"])+"\n")
    
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
        reactivity = np.array(data[name][0]).astype(float)
        structure = data[name][2]
        auc = AUC(structure,reactivity)
        auc_pred = AUC(structure,reactivity_pred)
        nan_mask = np.isnan(reactivity)
        reactivity = reactivity[~nan_mask]
        reactivity_pred = reactivity_pred[~nan_mask]
        corr,p = spearmanr(reactivity_pred,reactivity)
        rmse = RMSE(reactivity_pred,reactivity) 
        fout.write("\t".join([name,str(corr),str(p),str(auc),str(auc_pred),str(rmse)])+"\n")
    fout.close()



if __name__ == "__main__":
    main(args)

