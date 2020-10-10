#!/usr/bin/env python
from __future__ import print_function
import sys
import argparse
import re
import importlib

def main():
    parser = argparse.ArgumentParser(description='Filter Matrix by Expression or NaN values')
    parser.add_argument('--input', '-i',default="-",help="Input path, read from stdin by default.")
    parser.add_argument('--output','-o',default="-",help="Output path, print to stdout by default.")
    parser.add_argument('--name','-n',help="Sequence id prefix in the fastq file")
    parser.add_argument('--convert','-c',action="store_true",default=False,help="Whether convert U to T")
    parser.add_argument('--indice','-r',action="store_true",default=False,help="Whether append the indice of a sequence to the sequence id")
    parser.add_argument('--extract','-e',help="How to extract name from name line in ct file, should be the path of a python file contain a callable")
    args = parser.parse_args()
    extract_name = False
    if args.extract is not None:
        extract_name = True
        extract = importlib.import_module(args.extract)
    if args.input == "-":
        fin = sys.stdin
    else:
        fin = open(args.input)
    if args.output == "-":
        fout = sys.stdout
    else:
        fout = open(args.output,"w")
    sequence = ""
    NR = 0
    redundant_line = 0
    for line in fin:
        line = line.strip()
        if line.startswith("#") or len(line)==0:
            continue
        if "Energy" not in line and "ENERGY" not in line:
            redundant_line += 1
            continue
        NR += 1
        fields = re.split("\s+",line)
        length = int(fields[0])
        #print("Length of the sequence is:{}".format(length),file=sys.stderr)
        name_preffix = "" if args.name is None else args.name
        name_indice = str(NR) if args.indice else ""
        name_suffix = extract.extract(line) if extract_name else ""
        name_fields = [name_preffix,name_indice,name_suffix]
        name_fields = [item for item in name_fields if len(item)>0]
        name =  ":-:".join(name_fields)
        sequence = ""
        for i in range(length):
            line = next(fin)
            fields = re.split("\s+",line.strip())
            sequence += fields[1]
        fout.write("> "+name+"\n")
        sequence = sequence.upper()
        if args.convert:
            sequence = sequence.replace("U","T")
        fout.write(sequence+"\n")
    
    if redundant_line>0:
        print(args.name,file=sys.stderr)




if __name__=='__main__':
    main()

