#!/usr/bin/env python
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(description='Convert stockholm format to dot bracket format.')
parser.add_argument('--input', '-i',required=True,help="Input path, read from stdin by default.")
parser.add_argument('--outdir','-o',required=True,help="Output path, print to stdout by default.")
args = parser.parse_args()


inSeed = args.input
outDir = args.outdir
with open(inSeed,encoding="ISO-8859-1") as f:
    content = None
    accession = 'NA'
    for lineno, line in tqdm(enumerate(f)):
        if line.startswith('# STOCKHOLM 1.0'):
            if content:
                content = ''.join(content)
                with open(outDir+"/{}.seed".format(accession),"w") as f:
                    f.write(content)
            content = []
        elif line.startswith('#=GF AC'):
            accession = line.split()[-1]
            
        content.append(line)

