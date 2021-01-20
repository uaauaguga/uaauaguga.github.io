from Bio import Entrez
import sys
import json
import time
import numpy as np
import argparse
import requests
from pyfaidx import Fasta
import os

def main():
    parser = argparse.ArgumentParser(description='Download genebank sequence')
    parser.add_argument('--query', '-q',help="Genebank id for download",required=True)
    parser.add_argument('--fasta','-o',help="Full length sequence in genebank",default="Rfam/seeds/CMfinder-set-full-length.fasta")
    parser.add_argument('--download','-d',action="store_true",help="whether perform downloading",default=False)
    args = parser.parse_args()

    Entrez.email = "jinyf16@mails.tsinghua.edu.cn"
    gbIds = np.unique(open(args.query).read().strip().split("\n"))
    print("Load downloaded sequences ...")
    fastaObj = Fasta(args.fasta)
    downloaded = list(fastaObj.keys())
    fastaObj.close()
    print("Done .")
    print("Total downloaded sequences: {}".format(len(fastaObj.keys())))
    print("{} query sequence".format(gbIds.shape[0]))
    gbIds = np.setdiff1d(gbIds,downloaded)
    print("{} query sequence are not downloaded".format(gbIds.shape[0]))
    N = 0
    for gbId in gbIds:
        if gbId.startswith("URS"):
            N += 1
    print("{} query sequence are in RNAcentral annotation".format(N))

    if not args.download:
        sys.exit(0)


    fout = open(args.fasta,"a")

    for gbId in gbIds:
        try:
            if gbId.startswith("URS"):
                continue
                RNAcentralId,taxo = gbId.strip().split("_")
                print("Start retriving {} from RNAcentral...".format(gbId),file=sys.stderr)
                content = requests.get("https://rnacentral.org/api/v1/rna/{}/{}".format(RNAcentralId,taxo),headers={"Accept":"application/json"}).text
                data = json.loads(content)
                entry = ">" + gbId + " " + data["description"]
                sequence = data["sequence"]
                print(entry,file=fout)
                print(entry)
                print(sequence,file=fout)
                print("Done.")
            else:
                print("Start retriving {} from ncbi nucleotide...".format(gbId),file=sys.stderr)
                handle = Entrez.efetch(db="nucleotide", id=gbId, rettype="fasta", retmode="text")
                content = handle.read().strip()
                contents = content.split("\n")
                entry = contents[0] + "\n" + "".join(contents[1:])
                print(entry,file=fout)
                print("Done.")
        except:
            print("Error retriving {}, skip ...".format(gbId),file=sys.stderr)
        time.sleep(0.5)
    fout.close()
    if os.path.exists(args.fasta+".fai"):
        os.remove(args.fasta+".fai")

if __name__ == "__main__":
    main()