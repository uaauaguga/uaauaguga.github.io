#!/usr/bin/env python


import pickle
import argparse
import re

import shaker.rna_tools.rna_io as rio
import shaker.simushape as sim
import shaker.rna_tools.util as util
from eden import graph as eg


parser = argparse.ArgumentParser(description='Simulate reactivity from dot bracket file')
parser.add_argument('--input', '-i',required=True,help="Input known structure in dot bracket notation")
parser.add_argument('--output','-o',required=True,help="Output reactivity path")
parser.add_argument('--model','-m',help="Trained model for simulation",default="data/reactivity/shaker-model.pkl")
args = parser.parse_args()


print("Load model ...")
with open(args.model, 'rb') as fmdl:
    model =  pickle.load(fmdl)
print("Done .")
fout = open(args.output,"w")

first_entry = True
with open(args.input) as fin:
    for line in fin:
        line = line.strip()
        if line.startswith(">"):
            name = line.replace(">","")
            print("Processing {} ...".format(name))
            line = next(fin)
            sequence = line.strip()
            line = next(fin)
            dbn = line.split(" ")[0].strip()
            graph = util.sequence_dotbracket_to_graph(sequence,dbn)
            embedding = eg.vertex_vectorize([graph])[0]
            reactivity = model.predict(embedding).reshape(-1)
            data = [name]+list(reactivity.astype(str))
            fout.write("\t".join(data)+"\n")
        else:
            continue
fout.close()

