#!/usr/bin/env python
import sys
from glob import glob 
from Bio import SeqIO
from Bio.Alphabet import generic_rna
import argparse

def main():
    parser = argparse.ArgumentParser(description='convert cmfinder motifs to bed format')
    parser.add_argument('--input-prefix','-i',type=str, required=True, help='input prefix of cmfinder motif')
    parser.add_argument('--output','-o',type=str, required=True, help="output bed file")
    args = parser.parse_args()
    fout = open(args.output,"w")
    records = []
    for stk_path in glob(args.input_prefix + ".motif.h*"):
        if ".temp" in stk_path:
            continue
        motif_id = stk_path[stk_path.rfind("motif")+len("motif")+1:]
        fstk = SeqIO.parse(stk_path, "stockholm",alphabet=generic_rna)
        for record in fstk:
            seq_id = record.id
            location, score = record.description.split("\t")
            start, end = location.split("..")
            start, end = int(start), int(end)
            assert start < end
            records.append((seq_id,start,end,motif_id,score))
    records = sorted(records, key = lambda x: (x[0], x[1],x[3]))
    fout.write("track itemRgb=On\n")
    colors = [(76, 114, 176),  (221, 132, 82),  (85,168,104),
              (196, 78,  82),  (129, 114, 179), (147, 120,96),
              (218, 139, 195), (204, 185, 116), (100, 181, 205)]

    for seq_id, start, end, motif_id, score in records:
        if "." in motif_id:
            color = (140,140,140)
        else:
            #color = colors[hash(motif_id)%9]
            n_1, n_2 = motif_id[1:].split("_")
            n_1, n_2 = int(n_1), int(n_2)
            color = colors[(n_1*n_2-1)%9]
        fout.write(f'{seq_id}\t{start}\t{end}\t{motif_id}\t{score}\t+\t{start}\t{end}\t"{color[0]},{color[1]},{color[2]}"\n')
    fout.close()
if __name__ == "__main__":
    main()
