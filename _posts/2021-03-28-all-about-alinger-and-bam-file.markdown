---
layout: post
title:  "Aligners and Output Formatting"
date:   2021-03-28 11:18:15 +0800
usemathjax: true
categories: jekyll update
---

- This post tries to **record** what I've **heard of** about the following questions. 
  - Use which aligner, which parameter, to what reference sequence, to map what data
  - How to manipulate output of different aligners, or equivalently, how to manipulate bam files
- 学会调包
- Note this field is developing rapidly, thing wrote here may quickly out of date


## Famous tools for pairwise sequence comparison
- Short reads aligner
  - DNA-seq aligner
    - [bwa](http://bio-bwa.sourceforge.net/)
    - [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  - RNA-seq aligner
    - [STAR](https://github.com/alexdobin/STAR)
    - [hisat2](http://daehwankimlab.github.io/hisat2/)
- Long reads aligner
  - [minimap2](https://lh3.github.io/minimap2/minimap2.html)
- General propose aligner
  - blast
  - blat
- Whole genome aligner
  - [LASTZ](http://www.bx.psu.edu/~rsharris/lastz/)
  - [MUMmer](http://mummer.sourceforge.net/)



## Short reads mapping: reference and parameters

### Reads preprocessing
#### Adaptor trimming / quality trimming
- If you perform global alignment, then adapter trimming is required
#### Duplication removal / UMI handling
- You may remove PCR duplicates before alignment or after alignment
  - Pros: reduce computational burden in reads mapping.
  - Cons: less flexible
    - post alignment duplication removal allows "mark duplicates" instead of "remove duplicates"

### Short read aligners

- What short read aligners do is actually assign some annotations to each reads, that is which segment of the reads are assigned to which location of the genome in which form

- In different aligner, same thing may have different name


#### **Map to genome / map to transcriptome**
- Seems most project prefer map reads to genome, but some downstream analysis tools ([salmon](https://combine-lab.github.io/salmon/), [rsem](https://deweylab.github.io/RSEM/)) requires the reads to be aligned to transcriptome coordinate
- [STAR](https://github.com/alexdobin/STAR) could direct project genome aligned reads to transcriptome coordinate
  
#### **Which sequence to add in reference**
- Impact of haplotype and decoy sequences


#### **Gapped alignment / Ungapped alignment**

- Whether the aligner allow deletions and insertions 
- The CIGAR string in bam file contains `I` and `D` operation 
- Most aligners support gapped alignments, but note [bowtie1](http://bowtie-bio.sourceforge.net/manual.html) only support ungapped, so its output cannot be used for indel calling 

  
#### **DNA-seq alignment / RNA-seq alignment**

- The major difference between DNA-seq mapping and RNA-seq mapping is RNA-seq aligner perform spliced alignment
- For spliced aware aligner, the output bam file contains `N` in cigar operation
- Spliced-aware aligner usually accept genome annotation to define known splice junction
- The mapping quality calculation in DNA-seq aligner and RNA-seq aligner can be different, as RNA-seq is seldom used for variant calling

  

#### **Local alignment / End to end alignment**
- The scoring is different
- End to end alignment requires both end of the reads to align to the genome, while local alignment only requires a substring of reads to align to genome
- Local alignment is more sensitive, but may lead to false hit
- Local alignment produce so called soft-clipped segments if ends of the reads is not aligned, that is `S` cigar operation
- If you don't trim adaptor, you have to use local alignment, but for short reads (small RNA sequencing or riboseq for example) it can be dangerous
  
  
#### **Multimapper**
- Secondary alignment
- The problem of biased assignment
  
#### **Chimeric alignment**
- Split reads
- Supplementary alignment
  - Note supplementary alignment and secondary alignment are different things

#### **Unmapped reads**

#### **Paired end reads specific situation**
- Discordant reads
  - See some discussion here <https://www.biostars.org/p/278412/>
- Singleton
- Unmapped reads



#### Some parameter specification for reference

- ENCODE RNA-seq
- TCGA RNA-seq
- TCGA DNA-seq
- arriba's recommend parameter for fusion gene detection

## Alignment formats and APIs

### Manipulate bam file
- Format specification of bam file
- Flag in bam file
- Unmapped reads

#### Duplication handling

- For different task, duplication should be ideally handled in different ways
  - UMI is always favorable
  - Single end data and paired end data
  - In RNA-seq for gene expression quantification, better not remove duplicate is UMI is not present, see [here]({% post_url 2021-03-31-Dive-into-RNA-seq-quantification %}) for detail
  - For CLIP-seq, remember to collapse duplicate before peak calling

- Tools
  - Picard Markduplicates
  - samtools rmdup
    - See following posts for differences between samtools rmdup and picard
      - http://broadinstitute.github.io/picard/faq.html
        - samtools rmdup do not remove interchromosomal duplicates
      - http://www.bio-info-trainee.com/2003.html
  - [samblaster](https://github.com/GregoryFaust/samblaster)
    - https://academic.oup.com/bioinformatics/article/30/17/2503/2748175
  - [sambamba](https://lomereiter.github.io/sambamba/index.html)


### Tools for handling alignment in other format



