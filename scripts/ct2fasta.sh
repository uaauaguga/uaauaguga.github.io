#!/bin/bash
#indir=known-structure/RNA-strand/ct
#outdir=known-structure/RNA-strand/fasta
indir=known-structure/bpRNA/ct
outdir=known-structure/bpRNA/fasta
i=0
for file in $(ls $indir);do
name=${file%.*} 
input=${indir}/${name}.ct
output=${outdir}/${name}.fa
cat ${input} | scripts/ct2fasta.py -n ${name} -e extract > ${output} 
if [ $(echo "$i%500" | bc ) -eq 0 ];then
echo -e "$i sequence processed ."
fi
i=$(echo "$i+1" | bc)
done

