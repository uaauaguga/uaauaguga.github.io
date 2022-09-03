#!/usr/bin/env python
import argparse
from Bio import AlignIO
def get_identity(seq_1, seq_2):
    assert len(seq_2) == len(seq_2)
    identical = 0
    L = len(seq_1)
    for a,b in zip(seq_1,seq_2):
        if a == b:
            if a == "-" or a == ".":
                L -= 1
            else:
                identical += 1
    return identical/L


def main():
    parser = argparse.ArgumentParser(description='Get sequence pairs with sequence similarity lower then a given threshold')
    parser.add_argument('--input', '-i',required=True,help="Input RNA alignment in stockhold format")
    parser.add_argument('--output','-o',required=True,help="Output path, print to stdout by default.")
    parser.add_argument('--threshold','-t',type=float,default=0.9,help="Only retain sequence pairs with similarity lower than this threshold")
    args = parser.parse_args()


    align = AlignIO.read(args.input, "stockholm")
    records = [(record.id,str(record.seq)) for record in align]
    n = len(records)

    lines = []

    for i in range(n):
        seq_id_1, seq_1 = records[i]
        for j in range(i+1,n):
            seq_id_2, seq_2 = records[j]
            similarity = get_identity(seq_1,seq_2)
            if similarity <= args.threshold:
                lines.append(f"{seq_id_1}\t{seq_id_2}\t{similarity}")

    if len(lines) > 0:
        fout = open(args.output,"w")
        for line in lines:
            print(line, file = fout)



if __name__ == "__main__":
    main()

