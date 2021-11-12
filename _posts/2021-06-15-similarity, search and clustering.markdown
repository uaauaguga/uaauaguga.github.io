---
layout: post
title:  "Efficient Searching and Clustering"
date:   2021-06-15 20:41:20 +0800
usemathjax: true
categories: jekyll update
---

## Content

- Distance measure

- Indexing and hashing for exact and approximate searching

- Clustering with pairwise distance

- Clustering without explicit pairwise distance calculation


## Distance measure

- metric 

- metric learning


## Indexing and hashing for exact and approximate searching

### Data Structure for Exact Search

- kd-tree
- hierarchical k-means tree


### [Locality sensitive hashing](https://en.wikipedia.org/wiki/Locality-sensitive_hashing)

- A good tutorial <http://infolab.stanford.edu/~ullman/mmds/ch3a.pdf>
- Also see [4 Pictures that Explain LSH - Locality Sensitive Hashing Tutorial](https://randorithms.com/2019/09/19/Visual-LSH.html)
- <http://www.stat.ucdavis.edu/~chohsieh/teaching/ECS289G_Fall2016/lecture10.pdf>

- <https://chanzuckerberg.github.io/ExpressionMatrix2/doc/LshSlides-Nov2017.pdf>


#### Approximation for Hamming Distance
- Bit sampling LSH, that is random choose a subset of bits

#### Approximation for Jaccard Distance

- **minhash** can be utilized to approximate **Jaccard distance** between k-mer sets of two large genome, see <https://uaauaguga.github.io/jekyll/update/2021/05/15/Alignment-free-method-for-sequence-comparison.html>. The same technique can be applied to estimate the similarity between tow document, or large sparse binary vector.
- 我们如果随机取k行出来，近似的不是Jaccard距离，而是Hamming距离
- "min" here means set with smallest hash value
- If we encode k-mer sets with binary table, a minhash function h that maps sets to rows, or a implicit reordering of the table rows
- minhash直接输出的结果称为signature matrix,每一列是一个sample，每一行是不同的哈希值

#### Approximation for Cosine Distance

- "Signed Random Projections"


#### Reading
- 2014,[Hashing for Similarity Search: A Survey](https://arxiv.org/pdf/1408.2927.pdf)
- 2017,[Exact clustering in linear time](https://arxiv.org/pdf/1702.05425.pdf)

- <https://www.csc2.ncsu.edu/faculty/nfsamato/practical-graph-mining-with-R/slides/pdf/Graph_Cluster_Analysis.pdf>


## Clustering with pairwise distance

## Clustering without explicit pairwise distance calculation

### Repeatedly take advantage of locality sensitive hashing 
- <https://chanzuckerberg.github.io/ExpressionMatrix2/doc/LshSlides-Nov2017.pdf>

### Start from an affinity graph

- <https://web.iitd.ac.in/~bspanda/graphclustering.pdf>
- <https://www.csc2.ncsu.edu/faculty/nfsamato/practical-graph-mining-with-R/slides/pdf/Graph_Cluster_Analysis.pdf>
