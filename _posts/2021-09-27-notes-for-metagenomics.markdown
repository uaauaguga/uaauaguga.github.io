---
layout: post
title:  "Notes on Metagenomics"
date:   2021-09-27 20:14:51 +0800
usemathjax: true
categories: jekyll update
---

### Metagenomic profiling

- Assign reads to reference sequences
- Estimate abundance of different microbes from amplicon sequencing or shutgun sequencing

- alignment based methods vs. kmer based methods
  - [metaphlan2](https://huttenhower.sph.harvard.edu/metaphlan2/)
  - kraken2

- marker gene based methods vs. marker gene independent methods
  - metaphylan2, krakenunique
  - kraken2

- nucleotide level analysis vs. protein level analysis
  - kraken2
  - Kaiju

- benchmark study
  

### Assembly

- metaspades: memory intensive
  - <https://cab.spbu.ru/software/meta-spades/>
  - 2017, *Genome Research*,[metaSPAdes: a new versatile metagenomic assembler](https://genome.cshlp.org/content/27/5/824.long)

- megahit
  - <https://github.com/voutcn/megahit>
  - 2015, *Bioinformatics*, [MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph](https://academic.oup.com/bioinformatics/article/31/10/1674/177884)

- idba-ud
  - <https://github.com/loneknightpy/idba>
  - 2012, *Bioinformatics*, [IDBA-UD: a de novo assembler for single-cell and metagenomic sequencing data with highly uneven depth](https://academic.oup.com/bioinformatics/article/28/11/1420/266973)



### Strategy for removing contaminations

- 2014, *BMC Biology*, [Reagent and laboratory contamination can critically impact sequence-based microbiome analyses](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-014-0087-z)

- 2019, *Trends in Microbiology*, [Contamination in Low Microbial Biomass Microbiome Studies: Issues and Recommendations](https://www.sciencedirect.com/science/article/pii/S0966842X18302531)

- 2019, *Plos Computational Biology*, [Recentrifuge: Robust comparative analysis and contamination removal for metagenomics](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006967)


### Genome similarity

- mash distance
- whole genome alignment


### Distance metrics for microbe community
- [Bray-Curtis distance](https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity)
- [unifrac distance](https://en.wikipedia.org/wiki/UniFrac)
- weighted unifrac distance


### Community analysis based on pairwise distance

- [PCoA](https://en.wikipedia.org/wiki/Multidimensional_scaling) Analysis 
  - PCoA (principle coordinate analysis) and MDS (multi-dimesional scaling) are essentially the same thing (but this seems the term MDS is seldom used by the microbiology community ...)
  - project each sample in a distance matrix to a low dimensional space, attempt to preserve the pairwise distances
  - See 
    - <https://www.stat.pitt.edu/sungkyu/course/2221Fall13/lec8_mds_combined.pdf>
    - <https://repository.upenn.edu/cgi/viewcontent.cgi?article=1000&context=cis_papers>

- [permanova](https://en.wikipedia.org/wiki/Permutational_analysis_of_variance) analysis
  - For a given grouping (real grouping or shuffled grouping), a pseudo F statistics can be calculated
  - *P* value can be calculated by shuffling group labels 

### Differential analysis
- 2020, *Genome Biology*, [Assessment of statistical methods from single cell, bulk RNA-seq, and metagenomics applied to microbiome data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02104-1)

### Functional annotation
- eggnog: <http://eggnog.embl.de/>


### Visualization

- Krona plot
  - 2011, *BMC Bioinformatics* [Interactive metagenomic visualization in a Web browser](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-385)
  - <https://github.com/marbl/Krona/tree/master/KronaTools>
  - <https://github.com/khyox/recentrifuge>


- Visualization of phylogenetic tree and taxonomy tree
  - <https://github.com/joey711/phyloseq>
  - ggtree
  - <http://etetoolkit.org/>
    ```conda install -c conda-forge pyqt=4```
  - sankey diagram


### Genome annotation
- CDS prediction for prokaryotes
- CDS prediction for eukaryotes

- <https://github.com/nextgenusfs/funannotate>

### Useful Resources
- sequence id to taxo-id mapping
- <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/>
- <https://github.com/joey711/phyloseq>


### Antibiotic resistence gene finding


### Plasmid sequence prediction

- 2018, *NAR*, [PlasFlow: predicting plasmid sequences in metagenomic data using genome signatures](https://academic.oup.com/nar/article/46/6/e35/4807335)


### A collection of workflows
- 2021, *Plos Computational Biology*, [Metagenomics workflow for hybrid assembly, differential coverage binning, metatranscriptomics and pathway analysis (MUFFIN)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008716) also see <https://github.com/RVanDamme/MUFFIN>
- 2016, *Genome Biology*, [IMP: a pipeline for reproducible referenceindependent integrated metagenomic and metatranscriptomic analyses](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1116-8)

### Reviews for reference
- 2018, *Nature Review Microbiology*, [Best practices for analysing microbiomes](https://www.nature.com/articles/s41579-018-0029-9)

- 2020, *Genome Biology*, [The promise and challenge of cancer microbiome research](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02037-9)



- No update for a month ... Write something to demonstrate I still alive despite for its insiganificance...



