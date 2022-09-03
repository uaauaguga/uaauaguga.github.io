#!/usr/bin/env python
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(description='Split Rfam seed alignments, each RNA family as a stockholm format file.')
parser.add_argument('--input', '-i',required=True,help="Input path")
parser.add_argument('--outdir','-o',required=True,help="Output path")
args = parser.parse_args()


with open(args.input) as f:
    content = None
    accession = 'NA'
    for line in f:
        if line.startswith('INFERNAL1/a'):
            if content:
                content = ''.join(content)
                print(f"Saving {accession} ...")
                with open(args.outdir + "/{}.cm".format(accession),"w") as f:
                    f.write(content)
            content = []
        elif line.startswith('ACC'):
            accession = line.split()[-1]
            
        content.append(line)

