#!/bin/bash
for rfam in $(cat data/Rfam/info/rfam-id-st50-lt10.txt);do
echo $rfam
zcat data/Rfam/full-alignments/source/${rfam}.fa.gz | scripts/remove-duplicated-fasta-record.py > data/Rfam/full-alignments/fasta/${rfam}.fa
cmalign data/Rfam/cm-models/${rfam}.cm data/Rfam/full-alignments/fasta/${rfam}.fa > data/Rfam/full-alignments/stockholm/${rfam}.stk
done
