---
layout: post
title:  "Searching and Clustering"
date:   2021-06-15 20:41:20 +0800
usemathjax: true
categories: jekyll update
---


## Distance / Similarity Measure

## Indexing and hashing for exact / approximate searching

### Data Structure for Exact Search

- [kd-tree](https://en.wikipedia.org/wiki/K-d_tree): binary tree for efficient search Euclidean space 
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
- 2021, *Genome Research*, [SHARP: hyperfast and accurate processing of single-cell RNA-seq data via ensemble random projection](https://genome.cshlp.org/content/30/2/205.long)
- <https://chanzuckerberg.github.io/ExpressionMatrix2/doc/LshSlides-Nov2017.pdf>


## Clustering methods
- See <https://en.wikipedia.org/wiki/Cluster_analysis>

### K-means

- Overview
  - Random choose $$K$$ point as centroid
  - Assign each data point to nearest cetroid
  - Update centroid

- Time complexicity: $$O(nK)$$, where $$n$$ is number of sample points, $$K$$ is nnumber of predifined cluster
- Spce complexicity: $$O(n+K)$$

### [Hierarchical clustering](https://en.wikipedia.org/wiki/Hierarchical_clustering)

- Agglomerasive: start from 1 cluster which conatins all sample points
  - For each round
    - Input is some clusters with pairwise similarity / affinity matrix
    - Merge two closest clusters and update affinity matrix
  - **single linkage**
    - the proximity of two clusters is the maximum  of the proximity between any two points in the two clusters
  - **complete linkage**
    - the proximity of two clusters is the minimum of the proximity between any two points in the two clusters
  - **average linkage** alias [**UPGMA**](https://en.wikipedia.org/wiki/UPGMA) (Unweighted Pair Group Method using Arithmetic average)
    - the proximity of two clusters is average of pairwise proximities between sample in two clusters 
  - [**WPGMA**](https://en.wikipedia.org/wiki/WPGMA) 
    - weighted version of UPGMA
  - **centroid**
    - the proximity of two clusters is proximity between the centroids of two clusters
  - **ward**
    - the proximity is the increase in the squared error that results when two clusters are merged
  - hierarchical agglomerative clustering is generally computationally expensive
    - for general cases, $$O(n^3)$$ in time, $$O(n^2)$$ in space

- Divisive: start from $$n$$ clusters, each contains a sample sample point


### [DBSCAN](https://en.wikipedia.org/wiki/DBSCAN)
- Automatically determinesthe number of cluster
- Determinsitic
- Time complexicity: $$O(n)$$ * time for query points with in a given radius



### Types of similarity graph

- Dense similarity graph

- $$\epsilon$$-neighborhood graph

- k-nearest neighbor graph

- mutual k-nearest neighbor graph

- Shared Nearest Neighbor (SNN) similarity graph

- See <https://people.csail.mit.edu/dsontag/courses/ml14/notes/Luxburg07_tutorial_spectral_clustering.pdf>

### Graph base clustering

- The relationship between clustering and community detection

- Graph-based clustering, for instance,maps the task of finding clusters to the task of partitioning a proximity graph into connected components

- <https://web.iitd.ac.in/~bspanda/graphclustering.pdf>
- <https://www.csc2.ncsu.edu/faculty/nfsamato/practical-graph-mining-with-R/slides/pdf/Graph_Cluster_Analysis.pdf>

- Bron–Kerbosch algorithm 
- leiden algorithm


## Performance evaluation




### Clustering in single cell data analysis

- 2019, *Nature Review Genetics*, [Challenges in unsupervised clustering of single-cell RNA-seq data](https://www.nature.com/articles/s41576-018-0088-9)




## The relationship between clustering and searching

- Finding nearest neighbors can require computing the pairwise distance between all points. Often clusters and their cluster prototypes can be found much more efficiently.




