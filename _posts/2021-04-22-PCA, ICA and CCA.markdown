---
layout: post
title:  "PCA, CCA and ICA"
date:   2021-04-17 10:37:52 +0800
usemathjax: true
categories: jekyll update
---

- Three related dimensional reduction technique with application in gene expression data analysis
- PCA: principle component analysis
- CCA: canonical correlation analysis
- ICA: individual components analysis

- Consider we have a gene expression matrix, rows are genes, columns are samples

## PCA
- A linear projection that maximize variations or minimize reconstruction loss
- Suppose we have $$N$$ samples, $$X_1$$,$$X_2$$,...,$$X_n$$, each samples is a $$d$$ dimensional vector
- Sample mean is:  $$\bar{X}=\frac{\sum_{1}^{n}X_{i}}{N}$$
- Given a projection vector $$u_1$$ with dimension $$d$$, the projected vector is $$\hat{X_i} = u_{1}^TX_i$$ 
- Mean of projected vector is  $$\frac{\sum_{1}^{n}\hat{X_{i}}}{N} = \frac{\sum_{1}^{n}u_{1}^TX_i}{N} = u_{1}^T\bar{X}$$

#### The projected variance maximization formulation
- Variance of projected vector is:

$$
\begin{align*}
  &\frac{1}{N}\sum_{1}^{n}(\hat{X_{i}}-\bar{\hat{X_{i}}})^2 = \frac{1}{N} \sum_{1}^{n}(u_{1}^TX_i - u_{1}^T\bar{X})^2 \\
= & \frac{1}{N} \sum_{1}^{n}u_{1}^T(X_i - \bar{X})u_{1}^T(X_i - \bar{X}) \\
= & \frac{1}{N} \sum_{1}^{n}u_{1}^T(X_i - \bar{X})(X_i - \bar{X})^Tu_{1} \\
= & u_{1}^T[\frac{1}{N} \sum_{1}^{n}(X_i - \bar{X})(X_i - \bar{X})^T] u_{1} \\
= & u_{1}^TSu_{1}
\end{align*}
$$

- Where $$S=\frac{1}{N} \sum_{1}^{n}(X_i - \bar{X})(X_i - \bar{X})^T$$ is the covariance matrix

- Here we want to maximize $$u_{1}^TSu_{1}$$

- We add a restrain $$u_{1}^Tu_{1}=1$$  

- We have the language multiplier:

 $$u_{1}^TSu_{1} + \lambda(1-u_{1}^Tu_{1})$$

 - Take the derivitive, we have

 $$2Su_{1} - 2\lambda u_{1} = 0$$

 $$Su_{1} = \lambda u_{1} $$

 $$u_{1}^TSu_{1} = \lambda u_{1}^Tu_{1} = \lambda$$

 - Hence the projection that maximize variance of projected vector is the eigen vector corresponds to the largest eigen value 

 - As $$S^H=S$$, $$S$$ is a [Hermitian matrix ](https://en.wikipedia.org/wiki/Hermitian_matrix)(埃米尔特矩阵)
 
 - Hermitian matrix has orthogonal eigenvectors for distinct eigenvalues

 - Other eigen vector with smaller eigen values are other projections orthogonal to $$u_1$$

#### Reference
- [Principal component analysis: a review and recent developments](https://royalsocietypublishing.org/doi/10.1098/rsta.2015.0202)

## CCA
- Given two random vectors $$X$$ and $$Y$$, CCA seeks for two linear combination of these vectors, to maximize their correlation

- Some tutorial
  - <http://graphics.stanford.edu/courses/cs233-20-spring/ReferencedPapers/CCA_Weenik.pdf>


- Applications in cross studies and cross species integration of single cell transcriptomic profile

- Some tools
  - [Cross Decomposition in sklearn](https://scikit-learn.org/stable/modules/cross_decomposition.html)
  - [R package CCA](https://cran.r-project.org/web/packages/CCA/index.html)



## ICA
- [Independent Component Analysis:Algorithms and Applications](https://www.cs.helsinki.fi/u/ahyvarin/papers/NN00new.pdf)




