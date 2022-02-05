---
layout: post
title:  "Alignment Free Sequence Comparison"
date:   2021-05-15 12:07:32 +0800
usemathjax: true
categories: jekyll update
---

## Alignment free method for large scale genomic sequence analysis



### Approximate jaccard index calculation
- minhash
  - A instance for locality sensitive hashing that approximately preserve jaccard distance. (There are also some LSH that preserve hamming distance, etc)
  - k-minimum values (KMVs) sketching is a widely used variant of minhash
    - We have two genome
    - We have a hash function
    - For each genome, we calculate the hash value of every q-gram, take k smallest hash values
    - The overlap bewteen two set of hash value approximate jaccard distance beween all k-mers in two genomes 
  - Useful tools
    - [mash](https://github.com/marbl/Mash)
    - [sourmash](https://github.com/dib-lab/sourmash)
    - The `sketch.sh` utility in [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/), also see this post <https://www.biostars.org/p/234837/>

- Some applications
  - 2015, *NBT*, [Assembling large genomes with single-moleculesequencing and locality-sensitive hashing](https://www.nature.com/articles/nbt.3238) Use min-hash to define anchor between noisy long reads for sequence assembly
  - 2020,*Genome Biology*,[Metalign: efficient alignment-based metagenomic profiling via containment min hash](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02159-0) For metagenomic profiling, use min-hash to reduce the size of database to perform sequence alignment with minimap2


### k-mer sampling

- Counting / indexing  all k-mers in a genome requires large memory
- We may sample some representative k-mers

- minimizer / winnowing
  - sample 1 k-mer in w consecutive k-mers
  - For two substring of length w-k+1, if their sequence is identical, the same k-mer should be sampled
  - So in addition to specifying w and k, priority of all possible k-mer should be specified, and k-mer with hihgest priority is took
  - This strategy for k-mer sampling is originially proposed in 2004, *Bioinformatics*, [Reducing storage requirements for biological sequence comparison](https://academic.oup.com/bioinformatics/article/20/18/3363/202143) k-mer sampled by this strategy is called **minimizer**.
  - k-mer sampling has very wide applications in sequence comparisons
    - minimap, kraken, ...
  - further reading
    - <https://www.cs.cmu.edu/~gmarcais/slides/PAG2019.pdf>
    - <https://homolog.us/blogs/bioinfo/2017/10/25/intro-minimizer/>
    - <https://homolog.us/blogs/bioinfo/2014/04/06/2014-the-year-of-minimizers-in-bioinformatics/>
    - <https://genomeinformatics.github.io/mashmap/>
    - <https://academic.oup.com/bioinformatics/article/37/17/2563/6162158?login=true>
    - 2017, Plos Computational Biology, [Designing small universal k-mer hitting sets for improved analysis of high-throughput sequencing](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005777)

- Select m k-mers with smallest hash values 
  - Same as using m k-mers that corresponds to minhash sketches
  - See 2018, *Nature Communication*, [Clustering huge protein sequence sets in linear time](https://www.nature.com/articles/s41467-018-04964-5)
  - extract kmers such that the same k-mers tend to be extracted from homologous sequences.
  - avoid positional clustering of selected k-mers in order to be sensitive to detect local homologies in every region of a sequence
  - we would like to extract k-mers that tend to be conserved between homologous sequences
  - A MinHash sketch of size 1 is equivalent to minimizer, see the [mash](https://github.com/marbl/Mash) paper <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x> 

- 不同长度序列的minhash signature的大小是固定的，而minimizer的数量和序列的长度是正相关的


### LSH in other fields
- In addition to estimate Jaccard distance with minhash, we could use other hash function or other distance estimation
- See <http://infolab.stanford.edu/~ullman/mmds/ch3a.pdf> 



### Approximate member query and efficient counting
- [Bloom filter](https://en.wikipedia.org/wiki/Bloom_filters_in_bioinformatics)
  - Query whether an element exists in a large set. May generate false positive, but never false negative.
  - An alternative to the memory intensive hash table
  - A m bits array, K hash functions
  - For a new instance, each hash function map input to one of the positions in 1..m 
- Some variant
  - counting bloom filter
  - spectral bloom filter
- [Quotient filter](https://en.wikipedia.org/wiki/Quotient_filter)


## Reading List
- 2019, *Annual Review of Biomedical Data Science*, [Sketching and Sublinear Data Structures in Genomics](https://www.annualreviews.org/doi/abs/10.1146/annurev-biodatasci-072018-021156)
- 2019, *Genome Biology*, [When the levee breaks: a practical guide to sketching algorithms for processing the flood of genomic data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1809-x)

