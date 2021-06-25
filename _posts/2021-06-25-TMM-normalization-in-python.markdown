---
layout: post
title:  "TMM normalization in python"
date:   2021-06-25 13:16:51 +0800
usemathjax: true
categories: jekyll update
---


- TMM method in [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) package is a popular method for RNA seq count matrix, normalization, here I implement a python version


### Brief description 
- consider n gene, m sample in a (n,m) matrix
- take upper 75% quantile of each sample in expression matrix, divide it by library size, called f75 (m,)
- take sample nearest to mean as reference, you can set any sample as reference
- calculate TMM by compare each sample to the reference sample
- define: gene-wise log-fold-changes, also called log ratio of expression

$$M_g = \log_2{\frac{observation}{observation.libsize}} - \log_2{\frac{reference}{reference.libsize}}$$

- define: absolute expression level

$$A_g = \frac{1}{2}(\log_{2}\frac{observation}{observation.libsize}+\log_2{\frac{reference}{reference.libsize}})$$

- define: estimated asymptotic variance

$$ \\ V_g  = \frac{observation.libsize-observation}{observation.libsize*observation}+\frac{reference.libsize-reference}{reference.libsize*reference} 
      \\ = \frac{1}{observation} - \frac{1}{observation.libsize} + \frac{1}{reference} - \frac{1}{reference.libsize}$$

- Keep expression in quantile  [0.3,0.7] and relative expression in quantile [0.3,0.7] 

- calculate average of $$M_g$$, weighted by $$1/V_g$$

$$f = \frac{\sum\frac{1}{V_g}M_g}{\sum\frac{1}{V_g}}$$

- The TMM factor is $$2^f$$

### A python implementation


```python
import pandas as pd
import numpy as np


# Get TMM factor of a sample, given the reference sample
def get_TMM_factor(count_observed,count_reference,logratioTrim=0.3,sumTrim=0.05,Acutoff=-1e10):
    assert count_observed.shape[0] == count_reference.shape[0]
    
    # Mask genes with zero count in at least one of two samples
    mask = (count_observed>0)&(count_reference>0)
    count_observed = count_observed[mask].copy()
    count_reference = count_reference[mask].copy()
    
    observed_size = count_observed.sum()
    reference_size = count_reference.sum()
    
    # Calculate log fold change
    M = np.log2(count_observed/observed_size) - np.log2(count_reference/reference_size)

    
    if(M.abs().max()) < 1e-6:
        return 1
    

    # Calculate absolute expression 

    A = (np.log2(count_observed/observed_size) + np.log2(count_reference/reference_size))/2
    V =  1/count_observed - 1/observed_size +  1/count_reference -1/reference_size

    mask = (M < np.inf) & (A<np.inf) & (A>-1e10)
    M, A = M[mask], A[mask]
    

    # Get genes that is believed to not be differentially expressed in two sample
    # And genes which expression is not too high or too low

    n = M.shape[0]
    M_l, A_l = int(logratioTrim*n), int(sumTrim*n) 
    M_h, A_h = n - M_l,n - A_l
    
    M_rank,A_rank = M.rank(),A.rank()
    
    genes_selected = M_rank[(M_rank>=M_l) & (M_rank<=M_h) & (A_rank>=A_l) & (A_rank<=A_h)].index
    
    M,V = M.loc[genes_selected],V.loc[genes_selected]
    
    f = np.dot(1/V.values,M.values)/(1/V.values).sum()

    return 2**f


# Get TMM scaling factor for a expression matrix
def get_scaling_factor(expression_matrix,reference):
    scaling_factors = {}
    for col in expression_matrix.columns:
        f = get_TMM_factor(expression_matrix[col],expression_matrix[reference])
        scaling_factors[col] = f
    scaling_factors = pd.Series(scaling_factors)
    scaling_factors = scaling_factors/np.exp(np.log(scaling_factors).mean())
    return scaling_factors.loc[expression_matrix.columns]



# example usage
def example():
    expression_matrix = pd.read_csv("TCGA-CRC-counts.txt",sep="\t",index_col=0)
    expression_matrix = expression_matrix[expression_matrix.sum(axis=1)>0]
    scaling_factor = get_scaling_factor(expression_matrix,"TCGA-A6-5662-11A")

```


### Reference 

- TMM's publication: <https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25>
- Source code in edgeR: <https://rdrr.io/bioc/edgeR/src/R/calcNormFactors.R>