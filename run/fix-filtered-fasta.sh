#!/bin/bash
for data in RNAstrand combined;do 
mkdir -p data/reference-structure/model-organism/${data}/fasta-filtered-fixed
for organism in arabidopsis-thaliana  danio-rerio  drosophila-melanogaster  e-coli  homo-sapiens  mus-musculus  oryza-sativa  saccharomyces-cerevisiae;do
path=data/reference-structure/model-organism/${data}/fasta-filtered/${organism}.fa
output=data/reference-structure/model-organism/${data}/fasta-filtered-fixed/${organism}.fa
if [ -f $path ];then
echo -e "${organism}-${data}"
cat $path | awk '$0~/^>/{print;next;}{gsub("[^ACGTU]","N",$0);print;}' > ${output}
fi
done
done
