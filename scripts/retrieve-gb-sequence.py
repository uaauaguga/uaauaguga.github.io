from Bio import Entrez
import sys
import time
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Download genebank sequence')
parser.add_argument('--query', '-q',help="Genebank id for download",required=True)
parser.add_argument('--output','-o',help="Full length sequence in genebank",required=True)
parser.add_argument('--downloaded','-d',help="genebank id of sequences already downloaded",default="/data/jinyunfan/projects/structure-motif/Rfam/seeds/genebank-downloaded-ids.txt")
args = parser.parse_args()

Entrez.email = "jinyf16@mails.tsinghua.edu.cn"
gbIds = np.unique(open(args.query).read().strip().split("\n"))
downloaded = np.unique(open(args.downloaded).read().strip().split("\n"))
print("{} query sequence".format(gbIds.shape[0]))
gbIds = np.setdiff1d(gbIds,downloaded)
print("{} query sequence are not downloaded".format(gbIds.shape[0]))

fout = open(args.output,"w")

for gbId in gbIds:
    print("Start retriving {} ...".format(gbId),file=sys.stderr)
    try:
        handle = Entrez.efetch(db="nucleotide", id=gbId, rettype="fasta", retmode="text")
        content = handle.read().strip()
        print(content,file=fout)
    except:
        print("Error retriving {}, skip ...".format(gbId),file=sys.stderr)
        time.sleep(1)
fout.close()
