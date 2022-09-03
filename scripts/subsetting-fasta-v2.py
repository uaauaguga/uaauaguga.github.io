#!/usr/bin/env python
import argparse
from pyfaidx import Fasta
from tqdm import tqdm
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('subsetting fasta')


def main():
    parser = argparse.ArgumentParser(description='get sequence from fasta file')
    parser.add_argument('--input','-i',type=str,required=True,help="input fasta sequence")
    parser.add_argument('--output','-o',type=str,required=True,help="output fasta sequence")
    parser.add_argument('--seq-ids','-s',type=str,required=True,help="sequence ids to select")
    args = parser.parse_args()

    seq_ids = open(args.seq_ids).read().strip().split("\n")
    seq_ids = set(seq_ids)

    logger.info("load fasta ...")
    fasta = Fasta(args.input)
    fout = open(args.output,"w")
    
    logger.info("get sequences ...")
    width = 70
    for seq_id in tqdm(seq_ids):
        assert seq_id in fasta
        sequence = str(fasta[seq_id])
        fout.write(f">{seq_id}\n")
        p = 0
        while True:
            fout.write(sequence[p:p+width]+"\n")
            p += width
            if p >= len(sequence):
                break
    fout.close()
    logger.info("all done .")

if __name__ == "__main__":
    main()
    
