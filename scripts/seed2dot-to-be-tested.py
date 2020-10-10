#!/usr/bin/env python
import argparse
import sys
import re
from tqdm import tqdm
import numpy as np

parser = argparse.ArgumentParser(description='Extract concensus structure from seed alignment')
parser.add_argument('--input', '-i',required=True,help="Input path, read from stdin by default.")
parser.add_argument('--outdir','-o',required=True,help="Output path, print to stdout by default.")
parser.add_argument('--collapse','-cl',action="store_true",default=False,help="Collapse structure annotation.")
parser.add_argument('--convert','-cv',action="store_true",default=False,help="Convert U to T.")
args = parser.parse_args()

infile =  open(args.input)


seqDict = {}
print("Load seeds ...")
for line in tqdm(infile):
    line = line.strip()
    if len(line)==0 or line.startswith("/"):
        continue
    if line.startswith("#"):
        match = re.match(r"#=GC SS_cons\s+(.+)$",line)
        if match is not None:
            ss = match.group(1)
        continue
    name,seq = re.split("\s+",line)
    seqDict[name] = seq
print("Done .")

if args.collapse:
    ss = re.sub("[<(\[{]","(",ss)
    ss = re.sub("[>)\]}]",")",ss)
    ss = re.sub("[,\-.:_A-Ga-g]",".",ss)
ss =  np.array(list(ss))


print("Write dot files ...")
for name,seq in tqdm(seqDict.items()):
    if args.convert:
        seq  = seq.replace("U","T")
    seq = np.array(list(seq))
    idx = seq!="-"
    outfile = open(args.outdir+"/"+name.replace("/",":")+".dot","w")
    print(">"+name,file=outfile)
    print("".join(seq[idx]),file=outfile)
    print("".join(ss[idx]),file=outfile)
    outfile.close()
