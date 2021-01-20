#!/usr/bin/env python
import argparse
import os
from Bio import AlignIO

def getLength(path):
    totalLength = 0
    n = 0
    align = AlignIO.read(path, "stockholm")
    for record in align:
        n += 1
        sequence = str(record.seq).replace("-","")
        totalLength += len(sequence)
    return int(totalLength/n),n

def main():
    parser = argparse.ArgumentParser(description='Get length of alignments in stockhold format')
    parser.add_argument('--indir', '-i',required=True,help="Input dir contains stk files")
    parser.add_argument('--output','-o',required=True,help="Output table of average length")
    args = parser.parse_args()
    fout = open(args.output,"w")
    fout.write("rfam-id\tlength\tnumber\n")
    for stk in os.listdir(args.indir):
        rfamId = stk.split(".")[0]
        path = os.path.join(args.indir,stk)
        l, n  = getLength(path)
        fout.write(rfamId + "\t" + str(l) + "\t" + str(n) + "\n")
    fout.close()
    
   
    
if __name__ == "__main__":
    main()
