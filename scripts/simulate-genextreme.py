#!/usr/bin/env python
import argparse
import numpy as np
from scipy.stats import genextreme
from collections import defaultdict
import pandas as pd
import pickle
import re
from utilities import loadRecords,annotatePairing,writeRecords


def simulate(pairing,modelDict):
    reactivity = np.zeros(len(pairing))
    pairing = "UU" + pairing + "UU"
    for idx in range(len(reactivity)):
        structure = pairing[idx:idx+5]
        reactivity[idx] = modelDict[structure].rvs()
        if reactivity[idx] > 100:
            print(structure,reactivity[idx])
    reactivity[reactivity<0] = 0
    idx95 = reactivity.argsort()[int(reactivity.shape[0]*0.95)]
    reactivity[reactivity > reactivity[idx95]] = reactivity[idx95]

    # Apply a weighted sliding window smoothing
    smoothed = np.convolve(reactivity,np.array([1/6,2/3,1/6]), mode='valid')
    # Recover the original length of the array
    reactivity = np.concatenate([np.array([reactivity[0]]),smoothed,np.array([reactivity[-1]])])

    return reactivity

def main():
    parser = argparse.ArgumentParser(description='Simulate reactivity from dot bracket file')
    parser.add_argument('--input', '-i',required=True,help="Input known structure in dot bracket notation")
    parser.add_argument('--output','-o',required=True,help="Output reactivity path")
    parser.add_argument('--model','-m',required=True,help="Trained model for simulation")
    args = parser.parse_args()

    print("Load model ...")
    with open(args.model, 'rb') as f:
        modelDict =  pickle.load(f)

    print("Load input data ...")
    dataDict = loadRecords(args.input,data_type ="sequence,structure")
    dataDict = annotatePairing(dataDict)
   
    print("Perform simulation ...") 
    for seq_id in dataDict.keys():
        pairing = dataDict[seq_id]["pairing"]
        dataDict[seq_id]["reactivity"] = simulate(pairing,modelDict)
    print("Done .")

    writeRecords(dataDict,args.output,"reactivity")

if __name__ == "__main__":
    main() 
