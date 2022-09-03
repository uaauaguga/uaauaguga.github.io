#!/usr/bin/env python
import argparse

def main():
    parser = argparse.ArgumentParser(description='get sequence from fasta file')
    parser.add_argument('--input','-i',type=str,required=True,help="input fasta sequence")
    parser.add_argument('--output','-o',type=str,required=True,help="output fasta sequence")
    parser.add_argument('--seq-ids','-s',type=str,required=True,help="sequence ids to select")
    args = parser.parse_args()

    seq_ids = open(args.seq_ids).read().strip().split("\n")
    seq_ids = set(seq_ids)
    
    fin = open(args.input)
    fout = open(args.output,"w")
    
    for line in fin:
        if line.startswith(">"):
            seq_id = line.strip().split()[0][1:]
            kept = seq_id in seq_ids
            if kept:
                fout.write(line)
        else:
            if kept:
                fout.write(line)
    fin.close()
    fout.close()

if __name__ == "__main__":
    main()
    
