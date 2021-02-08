#!/bin/bash
indir=data/Rfam/seed-alignments
outdir=data/reference-structure/seed-alignments-structure
for file in $(ls ${indir});do
rfamid=${file%%.*}
echo $rfamid
scripts/stk2dbn.py -i ${indir}/${file} -o ${outdir}/${rfamid}.dot > log/get-known-structure/${rfamid}.log 2>&1
done
