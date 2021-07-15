---
layout: post
title:  "Split Sequencing Reads"
date:   2021-07-10 09:52:24 +0800
usemathjax: true
categories: jekyll update
---

- We may want to filter some unwanted sequences from fastq file prior to down stream analysis
  - rRNA for RNA seq data
  - host reads for metagenomic data 


## Remove rRNA sequence

- For single species RNA-seq, simply map reads to rRNA reference of these species, keep unaligned reads

- For metagenomic data, may try [sortmerna](https://github.com/biocore/sortmerna)

```bash
# download rRNA reference from their github repo https://github.com/biocore/sortmerna/tree/master/data/rRNA_databases
# build sortmerna index
sortmerna -ref reference/rRNA/fasta/rRNA.fasta -index 1 --idx-dir reference/rRNA/sortmerna-index 
# -index 1: build index and exit

# run sortmerna
 sortmerna -ref reference/rRNA/fasta/rRNA.fasta --idx-dir reference/rRNA/sortmerna-index --reads {input.fastq_1} --reads {input.fastq_2} --fastx --workdir {working_dir} --index 0 --out2 --other --threads {threads}
# -index 0: perform rRNA matching on existed index
```

- For annotate rRNA sequence in assembled sequence contigs of bacteria genome , try [barrnap](https://github.com/tseemann/barrnap)


## Remove host genome sequence

- align reads to host genome

- bmtagger
  - Used in HMP project, see <https://www.hmpdacc.org/hmp/doc/HumanSequenceRemoval_SOP.pdf>

- bbduk.sh in bbmap suites


