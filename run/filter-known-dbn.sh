#!/bin/bash
data=RNAstrand #Rfam
mkdir -p data/reference-structure/model-organism/${data}/dbn-filtered
for path in $(ls data/reference-structure/model-organism/${data}/fasta-filtered/*fa);do
file=${path##*/}
organism=${file%.*}
echo $organism
dbn=data/reference-structure/model-organism/${data}/dbn/${organism}.dot
echo $path
scripts/get-dbn-records.py -i ${dbn} -q $(cat $path | grep  '^>' | sed 's/^>//' | tr '\n' ',') -o data/reference-structure/model-organism/${data}/dbn-filtered/${organism}.dot
done
