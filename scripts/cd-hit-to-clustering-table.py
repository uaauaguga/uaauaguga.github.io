#!/usr/bin/env python
import argparse
import re
from tqdm import tqdm

def main():
    parser = argparse.ArgumentParser(description='Reformat cd-hit clustering table')
    parser.add_argument('--input', '-i', required=True, help="Input path")
    parser.add_argument('--output','-o', required=True, help="Output path")
    args = parser.parse_args()
    fin = open(args.input)
    fout = open(args.output,"w")
    # 0       190nt, >32711|RF00174-Cobalamin::OTU-16336:3300026516.a:Ga0209573_1000185:139054-139618(+)/219-408... *
    # 1       114nt, >RF01482-AdoCbl_riboswitch::OTU-12326:3300021604.a:Ga0226835_1000036:109575-110162(+)/271-384... at +/83.33% 
    pattern = r"(\d+)nt, >(.+)\.\.\. (.+)" 
    seq_ids = []
    centroid_id = ""
    for line in tqdm(fin):
        if line.startswith(">"):
            if len(seq_ids) > 0:
                assert len(centroid_id) > 0
                for seq_id in seq_ids:
                    fout.write(f"{seq_id}\t{centroid_id}\n")
            seq_ids = []
            centroid_id = "" 
        else:
            cluster_id, description  = line.strip().split("\t")
            m = re.match(pattern, description) 
            length, seq_id, centroid = m.groups()
            seq_ids.append(seq_id)
            if centroid == "*":
                centroid_id = seq_id
    if len(seq_ids) > 0:
        assert len(centroid_id) > 0
        for seq_id in seq_ids:
            fout.write(f"{seq_id}\t{centroid_id}\n")
    fin.close()
    fout.close()
if __name__ == "__main__":
    main()    
