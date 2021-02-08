#!/usr/bin/env python
import gzip
from tqdm import tqdm

seq_ids = open("data/accesion2taxid/seed-seq-ids.txt").read().strip().split("\n")
seq_ids = set([seq_id.split(".")[0] for seq_id in seq_ids])


path = "data/accesion2taxid/nucl_gb.accession2taxid.gz"
used = "data/accesion2taxid/seed-id2taxo.txt"

fout = open(used,"w")



with gzip.open(path) as f:
    _ = next(f)
    for line in tqdm(f):
        line = line.decode().strip()
        fields = line.split("\t")
        accession = fields[0]
        if accession in seq_ids:
            tax_id = fields[2]
            fout.write("{}\t{}\n".format(accession,tax_id))
            
fout.close()        
