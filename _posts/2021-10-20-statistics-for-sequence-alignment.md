---
layout: post
title:  "Statistics for Sequence Alignment"
date:   2021-10-20 14:34:51 +0800
usemathjax: true
categories: jekyll update
---

### Karlin-Altschul statistics
- Evaluate statistical significance of hits in blast search
- 1990, *PNAS*, S Karlin & S F Altschul, [Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes](https://www.pnas.org/content/87/6/2264.long)
- The expected number of **ungapped alignments** with score $$S$$ found with random sequences $$E$$ was demostrated to be
 
 $$E(S) = Kmne^{-\lambda S}$$

- $$m$$ is length of sequence 1, $$n$$ is length of sequence 2, $$K$$ is a constant

- Given the substitution matrix $$S_{i,j}$$, and frequency of the character $$a_{i}$$ and $$a_{j}$$

$$\sum_{i,j=1}^{r}p_{i}p_{j}e^{s_{ij} \lambda}=1$$

- Alternatively, the alignment score $$S$$ follows [extreme value distribution](https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution) (EVD), where $$u=\frac{lnKmn}{\lambda}$$

$$P(S<x)=e^{-e^{-\lambda(x-u)}}=e^{-Kmne^{-{\lambda} x}}$$

$$P(S{\geq}x)=1-P(S<x)=1-e^{-Kmne^{-{\lambda} x}}=1-e^{-E(x)}$$

- $$E(S)$$ can be interpreted as expected value in poisson distribution


- Readings:
  - <http://www.bioinfo.org.cn/lectures/index-47.html>
  - <https://personal.utdallas.edu/~prr105020/biol6385/2018/lecture/Stat_sig.pdf>
  - <http://pedagogix-tagc.univ-mrs.fr/courses/bioinfo_intro_prev/articles/sequence_alignment/Korf_BLAST_essential_OReilly.pdf>
  - <https://www.sciencedirect.com/science/article/pii/S0076687996660297?via%3Dihub>
