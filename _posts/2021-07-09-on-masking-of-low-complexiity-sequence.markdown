---
layout: post
title:  "On Masking of Low Complexity Sequence"
date:   2021-07-09 09:52:24 +0800
usemathjax: true
categories: jekyll update
---

### Handle low complexity sequence
- masking low complexity regions in a long sequence is sometimes favorable for database searching and motif finding, see <http://web.mit.edu/meme_v4.11.4/share/doc/glam2_tut.html>
- filter read with low sequence complexity in NGS data

### Method
- entropy
  - The [SEG algorithm](https://kodomo.fbb.msu.ru/FBB/year_10/ppt/SEG-93.pdf)
  - entropy of 3-mer in sliding window
- dust score of 3 mer in sliding window
  - see <https://kodomo.fbb.msu.ru/FBB/year_09/ppt/DUST.pdf>
- perform local alignment to sequence itself and known repeative elements
  

### Tools

- [repeatmasker](https://www.repeatmasker.org/), screens DNA sequences for interspersed repeats and low complexity DNA sequences

- NCBI's blast is shipped with multiple masking tools, see <https://www.ncbi.nlm.nih.gov/books/NBK569845/> for detail
  - segmasker
  - dustmasker
  - windowmasker

```{bash}
# example usage of dustmasker
dustmasker -outfmt fasta -parse_seqids -in {input.fasta} -out {masked.fasta}
# low complexity sequence set set to lower case, if you want to set it to N, run
sed '/^>/! s/[[:lower:]]/N/g' {masked.fasta} > {hardmasked.fasta} # see https://www.biostars.org/p/13677/s
```

- [dust](https://meme-suite.org/meme/doc/dust.html), a program shipped with MEME suites

```{bash}
dust sequences.fasta {cutoff} > sequences.masked.fasta
```

- NCBI's [seg](ftp://ftp.ncbi.nih.gov/pub/seg/seg/), segment sequence by local complexity. A very old tool ...

- sequence complexity (entropy and dust) based read filter in [prinseq](http://prinseq.sourceforge.net/), also see <https://chipster.csc.fi/manual/prinseq-complexity-filter.html>

- I write a small script that assign a single 3-mer count based dust or entropy score to each reads in fastq file, see <https://github.com/uaauaguga/NGS-Analysis-Notes/blob/master/scripts/sequence-complexity.py>


- 2011, *NAR*, [A new repeat-masking method enables specific detection of homologous sequences](https://academic.oup.com/nar/article/39/4/e23/1006710)
