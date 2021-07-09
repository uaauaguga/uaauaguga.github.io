---
layout: post
title:  "On Masking of Sequence with Low Comlexity "
date:   2021-07-09 09:52:24 +0800
usemathjax: true
categories: jekyll update
---

### handle low complexity sequence in biological sequence analysis
- masking low complexity regions in a long sequence
- filter read with low sequence complexity in NGS data


### Repeats vs. Low complexity sequence



### Method
- entropy
  - The [SEG algorithm](https://kodomo.fbb.msu.ru/FBB/year_10/ppt/SEG-93.pdf)
  - entropy of 3-mer
- dust score of 3 mer
  - see <https://kodomo.fbb.msu.ru/FBB/year_09/ppt/DUST.pdf>
  - for long sequence, local complexity is determined by calculating dust / entropy in sliding window with given step length and window size 

### Tools

- segmasker, a repeat masking tool shipped with NCBI's blast

```bash
segmasker -in sequences.fasta -outfmt fasta > sequences.masked.fasta
```

- [dust](https://meme-suite.org/meme/doc/dust.html), a program shipped with MEME suites

```{bash}
dust sequences.fasta {cutoff} > sequences.masked.fasta
```

- NCBI's [seg](ftp://ftp.ncbi.nih.gov/pub/seg/seg/), segment sequence by local complexity. A very old tool ...

- sequence complexity (entropy and dust) based read filter in [prinseq](http://prinseq.sourceforge.net/), also see 

- I write a small script that assign a single dust or entropy score to each reads in fastq file
