#!/bin/bash
for data in  RNAstrand combined;do 
mkdir -p data/reference-structure/model-organism/${data}/bowtie2-index
for organism in arabidopsis-thaliana  danio-rerio  drosophila-melanogaster  e-coli  homo-sapiens  mus-musculus  oryza-sativa  saccharomyces-cerevisiae;do
path=data/reference-structure/model-organism/${data}/fasta-filtered-fixed/${organism}.fa
if [ -f $path ];then
echo -e "${data}-${organism}"
#rm -r data/reference-structure/model-organism/${data}/star-index/${organism}
#scripts/build-index.py --fasta ${path} --aligner star --prefix data/reference-structure/model-organism/${data}/star-index/${organism}/ --tmp-dir tmp/star/${data}-${organism} > log/build-star-index/${data}-${organism}.log 2>&1
scripts/build-index.py --fasta ${path} --aligner bowtie2 --prefix data/reference-structure/model-organism/${data}/bowtie2-index/${organism}/sequence  > log/build-bowtie2-index/${data}-${organism}.log 2>&1
fi
done
done
