---
layout: post
title:  "Multiple Test Correction"
date:   2021-03-20 1:11:20 +0800
usemathjax: true
categories: jekyll update
---

## Manipulate Genomic Features

- Genomic data is often genomic ranges associate with some metadata
- gtf, bed, vcf, sam, wig, bedgraph ...
- Multiple tools is available for manipulate such genomic ranges

### A closer look at file formats

- The coordinate system
- gtf file
- gff file
- bed and bed12 format

## Useful tools

- Command line tools
  - bedtools
- Bioconductor package
  - GenomicRanges

## Some task related to genome range operation

#### Take intersection of feature sets
- Annotate input features with known feature set
- Command line
- R
- Take some arithmetic operation across genomic interval
- Command line
- R

#### Take union

- Merge bed file
- Merge / reduce exon from different isoforms from a gene

#### Take difference

- Get intron location from gff/gtf file

### Algorithm behind these tools

- Interval tree

### Implement similar function in other programming language

- perl
- python
- C++
  - https://github.com/lh3/cgranges





