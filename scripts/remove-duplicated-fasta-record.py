#!/usr/bin/env python
import argparse
import sys
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='Remove sequence with duplicatedname')
    parser.add_argument('--input', '-i',help="Input fasta file",default="-")
    parser.add_argument('--output','-o',help="Output fasta file",default="-")
    parser.add_argument('--number','-n',type=int,help="Output fasta file",default=2000)
    args = parser.parse_args()
    fin = sys.stdin if args.input == "-" else open(args.input) 
    fout = sys.stdout if args.output == "-" else open(args.input,"w")
    sequences = defaultdict(str)
    for line in fin:
        line = line.strip()
        if line.startswith(">"):
            name = line.replace(">","")
        else:
            sequences[name] += line
        if len(sequences.keys()) > args.number:
            break
    for name,sequence in sequences.items():
        fout.write(">"+name+"\n")
        fout.write(sequence+"\n")
    fin.close()
    fout.close()

if __name__ == "__main__":
    main()
            
            
    

