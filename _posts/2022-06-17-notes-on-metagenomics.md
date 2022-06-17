---
layout: post
title:  "Sequence Motif Modeling Related"
date:   2022-06-17 15:15:41 +0800
usemathjax: true
categories: jekyll update
---

- Estimate abundance of different microbes from amplicon sequencing or shutgun sequencing

- alignment based methods vs. kmer based methods
  - metaphlan
  - kraken
  - [metaphlan2](https://huttenhower.sph.harvard.edu/metaphlan2/)
  - kraken2

- marker gene based methods vs. marker gene independent methods
  - metaphylan, krakenunique
  - kraken
  - metaphylan2, krakenunique
  - kraken2

- nucleotide level analysis vs. protein level analysis
  - kraken
  - kraken2
  - Kaiju

- benchmark study


### Assembly

- metaspades: memory intesive
@@ -43,6 +46,12 @@ categories: jekyll update

### Strategy for removing contaminations

- 2014, *BMC Biology*, [Reagent and laboratory contamination can critically impact sequence-based microbiome analyses](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-014-0087-z)

- 2019, *Trends in Microbiology*, [Contamination in Low Microbial Biomass Microbiome Studies: Issues and Recommendations](https://www.sciencedirect.com/science/article/pii/S0966842X18302531)

- 2019, *Plos Computational Biology*, [Recentrifuge: Robust comparative analysis and contamination removal for metagenomics](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006967)


### Genome similarity

@@ -55,25 +64,31 @@ categories: jekyll update
- [unifrac distance](https://en.wikipedia.org/wiki/UniFrac)
- weighted unifrac distance

- PCoA Analysis 

### Community analysis based on pairwise distance

- [PCoA](https://en.wikipedia.org/wiki/Multidimensional_scaling) Analysis 
  - PCoA (principle coordinate analysis) and MDS (multi-dimesional scaling) are essentially the same thing (but this seems the term MDS is seldom used by the microbiology community ...)
  - project each sample in a distance matrix to a low dimensional space, attempt to preserve the pairwise distances
  - See 
    - <https://en.wikipedia.org/wiki/Multidimensional_scaling>
    - <https://www.stat.pitt.edu/sungkyu/course/2221Fall13/lec8_mds_combined.pdf>
    - <https://repository.upenn.edu/cgi/viewcontent.cgi?article=1000&context=cis_papers>

- [permanova analysis](https://en.wikipedia.org/wiki/Permutational_analysis_of_variance)
- [permanova](https://en.wikipedia.org/wiki/Permutational_analysis_of_variance) analysis
  - For a given grouping (real grouping or shuffled grouping), a pseudo F statistics can be calculated
  - *P* value can be calculated by shuffling group labels 

### Differential analysis

- 2020, *Genome Biology*, [Assessment of statistical methods from single cell, bulk RNA-seq, and metagenomics applied to microbiome data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02104-1)

### Functional annotation
- eggnog: <http://eggnog.embl.de/>


### Visualization



### Genome annotation
- CDS prediction for prokaryotes
- CDS prediction for eukaryotes


## Useful resources

- <https://github.com/ngs-docs/2017-ucsc-metagenomics>
- <https://github.com/GenomicsAotearoa/metagenomics_summer_school>