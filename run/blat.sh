#!/bin/bash
for data in RNAstrand combined;do #Rfam
mkdir -p data/reference-structure/model-organism/${data}/blat
for organism in arabidopsis-thaliana  danio-rerio  drosophila-melanogaster  e-coli  homo-sapiens  mus-musculus  oryza-sativa  saccharomyces-cerevisiae;do
path=data/reference-structure/model-organism/${data}/fasta-filtered/${organism}.fa
if [ -f $path ];then
bsub -q TEST-1U "blat genome/${organism}/fasta/genome.fasta $path -q=dnax -t=dnax -stepSize=11 -ooc=genome/${organism}/ooc/genome-11.ooc data/reference-structure/model-organism/${data}/blat/${organism}.psl > log/blat/${organism}-${data}.log 2>&1"
fi
done
done
