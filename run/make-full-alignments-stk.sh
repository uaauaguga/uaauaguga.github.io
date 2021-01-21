#!/bin/bash
for rfam in $(cat data/Rfam/info/rfam-id-st50-lt10.txt);do
#gunzip data/Rfam/full-alignments/fasta/${rfam}.fa.gz 
echo $rfam
bsub -q TEST-1U "cmalign data/Rfam/cm-models/${rfam}.cm data/Rfam/full-alignments/fasta/${rfam}.fa > data/Rfam/full-alignments/stockholm/${rfam}.stk"
done
