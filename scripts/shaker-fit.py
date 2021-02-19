#!/usr/bin/env python
import pickle
import argparse

import shaker.rna_tools.rna_io as rio
import shaker.simushape as sim
import shaker.rna_tools.util as util
from eden import graph as eg


parser = argparse.ArgumentParser(description='Train SHAKER model')
parser.add_argument('--dot-bracket', '-dbn',required=True, help="Input dot bracket notation file")
parser.add_argument('--reactivity','-r',required=True, help="Path of the reactivity")
parser.add_argument('--model','-m',required=True,help="Path of the output model")
args = parser.parse_args()


def main(args):
    print("Load input data ...")
    data = rio.get_all_data(args.reactivity,args.dot_bracket) 
    print("Done .")
    print("Train SHAKER model ...")
    model  = sim.make_model(data,data.keys())
    print("Done .")
    print("Saving the model...")
    with open(args.model, 'wb') as f:
        pickle.dump(model, f)
    print("Done .")
    


if __name__ == "__main__":
    main(args)

