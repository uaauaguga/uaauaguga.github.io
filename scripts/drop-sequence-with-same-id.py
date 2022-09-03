#!/usr/bin/env python
import argparse

def main():
    parser = argparse.ArgumentParser(description='drop sequence with identical sequence')
    parser.add_argument('--input','-i',type=str, required=True, help='input fasta file')
    parser.add_argument('--output','-o',type=str, required=True, help="output fasta file")
    args = parser.parse_args()
    
    fout = open(args.output,"w")
    seq_ids = set()
    with open(args.input) as f:
        for line in f:
            if line.startswith(">"):
                fields = line.strip()[1:].split(" ")
                seq_id = fields[0]
                if seq_id in seq_ids:
                    skipped = True
                else:
                    skipped = False
                    seq_ids.add(seq_id)
            if not skipped:
                fout.write(line)
        



if __name__ == "__main__":
    main()
