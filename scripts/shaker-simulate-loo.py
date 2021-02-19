#!/usr/bin/env python
import pickle
import argparse

import shaker.rna_tools.rna_io as rio
import shaker.simushape as sim
import shaker.rna_tools.util as util
from eden import graph as eg
from scipy.stats import spearmanr
import numpy as np

parser = argparse.ArgumentParser(description='Train SHAKER model')
parser.add_argument('--dot-bracket', '-dbn',required=True, help="Input dot bracket notation file")
parser.add_argument('--reactivity','-r',required=True, help="Path of the reactivity")
parser.add_argument('--output','-o',required=True,help="Result of cross validation")
args = parser.parse_args()

def RMSE(predictions, targets):
    predictions = np.array(predictions)
    targets = np.array(targets)
    return np.sqrt(((predictions - targets) ** 2).mean())


def main(args):
    print "Load input data ..."
    data = rio.get_all_data(args.reactivity,args.dot_bracket) 
    print "Done ."
    print "Train SHAKER model ..."
    
    fout = open(args.output,"w")
    fout.write("\t".join(["name","spearmanr","p-value","RMSE"])+"\n")
    for name in data.keys():
        keys = set(data.keys())
        keys.remove(name)
        model  = sim.make_model(data,list(keys))
        graph = util.sequence_dotbracket_to_graph(data[name][1],data[name][2])
        embedding = eg.vertex_vectorize([graph])[0]
        reactivity_pred = model.predict(embedding).reshape(-1)
        reactivity = np.array(data[name][0]).astype(float)
        reactivity[np.isnan(reactivity)] = 0
        corr,p = spearmanr(reactivity_pred,reactivity)
        rmse = RMSE(reactivity_pred,reactivity) 
        fout.write("\t".join([name,str(corr),str(p),str(rmse)])+"\n")
    fout.close()



if __name__ == "__main__":
    main(args)

