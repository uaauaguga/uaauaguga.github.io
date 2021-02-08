#!/bin/bash
for organism in danio-rerio;do #arabidopsis-thaliana homo-sapiens  mus-musculus  oryza-sativa  saccharomyces-cerevisiae;do
bsub -q TEST-1U "blat genome/${organism}/fasta/genome.fasta data/reference-structure/model-organism/fasta-filtered/${organism}.fa -q=dnax -t=dnax -stepSize=11 -ooc=genome/${organism}/ooc/genome-11.ooc data/reference-structure/model-organism/blat/${organism}.psl > log/blat/${organism}.log 2>&1"
done
