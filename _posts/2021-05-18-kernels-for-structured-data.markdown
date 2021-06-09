---
layout: post
title:  "Kernels for Structured Data"
date:   2021-05-18 10:21:39 +0800
usemathjax: true
categories: jekyll update
---

- We can use a transformation $$\phi$$ to map data space $$R^d$$ to a higher dimensional feature space $$R^p$$
- Kernel is a function that implicitly and efficiently calculate inner product in feature space
  - Calculate the inner product without calculate $$\phi$$
  - That means we don't have to know the explicit analytical form of $$\phi$$
  - Kernel function can be extended to structured data, like string, tree, and graph
    

 $$\kappa(x_i,x_j) = <\phi(x_i),\phi(x_j)^T> $$


### String kernel
- 2002, *Journal of Machine Learning Research 2*, [Text Classification using String Kernels](https://www.jmlr.org/papers/volume2/lodhi02a/lodhi02a.pdf)
  - String subsequence kernel
    - feature transformation
    - kernel function

### Tree kernel


### Graph kernel
- The R-convolution [Convolution kernels on discrete structures](http://www0.cs.ucl.ac.uk/staff/m.pontil/reading/haussler.pdf)
- Decompose two object into two sets of components, perform pairwise comparison between these components, and summarize to a numeric value

- Graphlet
- Subtree Patterns
  - Weisfeiler-Lehman algorithm
  - Weisfeiler-Lehman 
  - 2010, Journal of Machine Learning Research, [Weisfeiler-Lehman graph kernels](https://people.mpi-inf.mpg.de/~mehlhorn/ftp/genWLpaper.pdf)
- Random Walks

### Implementations
- [graphkit-learn](https://graphkit-learn.readthedocs.io/en/master/)
- [GraphKernels](https://github.com/BorgwardtLab/GraphKernels)
- [GraKeL](https://ysig.github.io/GraKeL/0.1a8/documentation.html)
- [EDeN](https://github.com/fabriziocosta/EDeN)
