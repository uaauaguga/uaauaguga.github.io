#!/usr/bin/env python
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Get reference structure of a specified species')
    parser.add_argument('--dotdir', '-d',help="Dir contains reference dbn files",default="data/reference-structure/seed-alignments-structure")
    parser.add_argument('--idtotax',help="Sequence id to tax id mapping",default="data/reference-structure/seed-id2taxo.txt")
    parser.add_argument('--output','-o',help="Output dbn files",required=True)
    parser.add_argument('--taxid','-t',help="Taxid of specified species",required=True)
    args = parser.parse_args()
    print("Load seed ids of species {} ...".format(args.taxid))
    used_seq_ids = set()
    with open(args.idtotax) as f:
        for line in f:
            seqid,taxid = line.strip().split("\t")
            if taxid == args.taxid:
                used_seq_ids.add(seqid)
    print("{} sequence ids loaded".format(len(used_seq_ids)))
    

    fout = open(args.output,"w")

    n = 0
    sequences = set()
    for file in os.listdir(args.dotdir):
        rfam_id = file.split(".")[0]
        path = os.path.join(args.dotdir,file)
        with open(path) as f:
            for line in f:
               assert line.startswith(">") 
               name = line.strip().replace(">","")
               gbid = name.split(":")[0].split(".")[0]
               sequence = next(f).strip()
               structure = next(f).strip()
               if sequence in sequences:
                   continue
               else:
                   sequences.add(sequence)
               if structure.count("(") == 0:
                   continue
               if gbid in used_seq_ids:
                   n += 1
                   name = ">" + rfam_id + ":" + name 
                   print(name,file=fout)
                   print(sequence.replace("U","T"),file=fout)
                   print(structure,file=fout)
    print("Get {} sequences with known reference structure".format(n))

if __name__ == "__main__":
    main()
               
    

