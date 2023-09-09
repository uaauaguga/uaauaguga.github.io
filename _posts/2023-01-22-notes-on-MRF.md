---
layout: post
title:  "Notes on Markov Random Field"
date:   2023-01-22 10:37:51 +0800
usemathjax: true
categories: jekyll update
---


# Notes on Markov Random Field

## Definition and parameterization

- Markov random field is an undirectional probablistic graph model
- Consider random variables $$X=(X_1,X_2,...,X_n)$$
- These random variables have conditional dependency defined by a graph $$G$$
- The graph $$G$$ have $$m$$ cliques (fully connected subgraphs), denoted as $$D_1,...,D_m$$
- Then distribution of the random variable can be expressed as product of clique potentials: 

$$
\large{P(X)=\frac{\prod_{i=1}^{m}\phi_{i}(D_i)}{Z}}
$$

- $$\phi_{i}(D_i)$$ is potential of clique $$D_i$$
- $$Z$$ is called partition function

$$Z = \sum_{X}\prod_{i=1}^{m}\phi_{i}(D_i)$$

- Calculation of $$Z$$ is generally computationally intractable

- It can also expressed as an log-linear model

$$
\large{P(X)=\frac{1}{Z}\exp{[-\sum^{k}_{i}w_{i}f_{i}(D_{i})]}}
$$

- $$w_{i}$$ are weights
- $$f_{i}(D_{i})$$ is the feature of complete subgraph $$D_{i}$$
- $$Z$$ is the normalization constant, or partition function



## Image denoising 


## Image segmentation


## Ising / Potts model and protein / RNA contact prediction 

### Ising model

- Random variable $$X_i$$ takes value of either 1 or -1
- The energy of edge (i,j) is 

$$\epsilon_{i,j}(x_i,x_j)=w_{i,j}x_{i}x_{j}$$

- The energy of variable $$X_i$$ is  

$$u_{i}x_{i}$$

- The distribution takes the form

$$P(X_1,...,X_n) = \frac{1}{Z}\exp(-\sum_{i<j}w_{i,j}x_{i}x_{j}-\sum_{i}u_{i}x_{i})$$

### Potts model and protein / RNA structure modeling

- A protein sequence in multiple sequence alignment: $$(\sigma_1,\sigma_2,\sigma_3,...,\sigma_N)$$
- $$\sigma_1$$ takes value from the protein/RNA alphabet and gap "-"


- mfDCA (mean field direct coupling analysis)
- plmDCA  (pseudolikelihood maximization direct coupling analysis)

### Restricted bolzmann machine

## Conditional random field

- Direct modeling the conditional distribution $$P(y\|X)$$, that is given observations of a subset of variables $$X$$, modeling the distribution remaining variables $$Y$$ conditioned on $$X$$
- Combining the advantages of discriminative classification and graphical modeling
- Combining the ability to compactly model multivariate outputs y with the ability to leverage
a large number of input features x for prediction

$$
\large{P(Y|X)=\frac{\prod_{i=1}^{m}\phi_{i}(D_i)}{Z}}
$$


$$Z = \sum_{Y}\prod_{i=1}^{m}\phi_{i}(D_i)$$

- Note $$X$$ is not involved in calculation of partition function

- Logistic regression is a special case of conditional random field



{:refdef: style="text-align: center;" }
![image]({{site.baseurl}}/images/2023-01-23-logistic-regression.png){:width="30%"}
{:refdef}

$$
\large{
    P(y_i|X) = \frac{\exp{(\beta_{y_{i}0}x_0+\sum_{1}^{m}\beta_{y_{i}k}x_k)}}{\sum_{y_{i}}\exp{(\beta_{y_{i}0}x_0+\sum_{1}^{m}\beta_{y_{i}k}x_k)}}
    = \frac{1}{Z}\exp{-\sum_{0}^{k}w_{y_{i}k}f_{k}(y,X)}
    }
$$


- It takes a form same as log-linear model of markov random field


## Linear chain CRF


$$
\large{
p(y\|X) = \frac{1}{Z(X)}\exp{[\sum_{k=1}^{K}\lambda_{k}f_{k}(y_t,y_{t-1},x_t)]}
}
$$

$$
\large{
Z(X) = \sum_{y}\exp{[\sum_{k=1}^{K}\lambda_{k}f_{k}(y_t,y_{t-1},x_t)]}
}
$$


## Reference

- [An Introduction to Conditional Random Fields](https://homepages.inf.ed.ac.uk/csutton/publications/crftut-fnt.pdf)
