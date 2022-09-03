#!/usr/bin/env python
import argparse
from utilities import loadRecords 
import os

def main():
    parser = argparse.ArgumentParser(description='Split dot bracket files')
    parser.add_argument('--dot-bracket','-dbn',type=str,required=True,help="Input dot bracket file")
    parser.add_argument('--outdir','-o',type=str,required=True,help="Output dir of the result.")
    args = parser.parse_args()
    dataDict = loadRecords(args.dot_bracket,"sequence,structure")
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    for seq_id in dataDict.keys():
        path = os.path.join(args.outdir,seq_id+".dot")
        with open(path,"w") as f:
            f.write(">"+seq_id+"\n")
            f.write(dataDict[seq_id]["sequence"]+"\n")
            f.write(dataDict[seq_id]["structure"].split(" ")[0]+"\n")
        

if __name__ == "__main__":
    main()        
