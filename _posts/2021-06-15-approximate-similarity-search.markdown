---
layout: post
title:  "Hash and Approximate Similarity Search"
date:   2021-06-15 20:41:20 +0800
usemathjax: true
categories: jekyll update
---

## Locality sensitive hashing

- A good tutorial <http://infolab.stanford.edu/~ullman/mmds/ch3a.pdf>


### Approximation for Jaccard Distance

- **minhash** can be utilized to approximate **Jaccard distance** between k-mer sets of two large genome, see <https://uaauaguga.github.io/jekyll/update/2021/05/15/Alignment-free-method-for-sequence-comparison.html>. The same technique can be applied to estimate the similarity between tow document, or large sparse binary vector.
- 我们如果随机取k行出来，近似的不是Jaccard距离，而是Hamming距离
- "min" here means set with smallest hash value
- If we encode k-mer sets with binary table, a minhash function h that maps sets to rows, or a implicit reordering of the table rows
- minhash直接输出的结果称为signature matrix,每一列是一个sample，每一行是不同的哈希值


### Approximation for Cosine Distance

- sketch 


### Reading
- 2014,[Hashing for Similarity Search: A Survey](https://arxiv.org/pdf/1408.2927.pdf)