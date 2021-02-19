#!/usr/bin/env python
import argparse

def main():
    parser = argparse.ArgumentParser(description='Get queried record from dot bracket file')
    parser.add_argument('--input', '-i',required=True,help="Input dot file")
    parser.add_argument('--query', '-q',required=True,help="Queried sequence ids, seperated by coma")
    parser.add_argument('--output','-o',required=True,help="Output dot file")
    args = parser.parse_args() 
    seq_ids = set(args.query.split(","))
    fin  = open(args.input)
    fout = open(args.output,"w")
    for line in fin:
        name = line
        seq_id = line.strip()[1:].split()[0].upper()
        sequence = next(fin)
        structure = next(fin)
        if seq_id in seq_ids:
            fout.write(name)
            fout.write(sequence)
            fout.write(structure)
    fin.close()
    fout.close()


if __name__ == "__main__":
    main() 
            
             
