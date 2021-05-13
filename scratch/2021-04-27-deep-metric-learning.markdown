---
layout: post
title:  "Deep Metric Learning"
date:   2021-04-30 19:38:33 +0800
usemathjax: true
categories: jekyll update
---

## Metric Learning

- A Metric is a function that quantifies a distance between every pair of elements in a set, thus inducing a measure of similarity.
- Metric fucntion can be defined, or learnt from data
  - Euclidian distance, Hamming distance, ...
  - [Mahalanobis Distance](https://en.wikipedia.org/wiki/Mahalanobis_distance) ($$\Sigma$$ is covariance matrix)
    - Same as whitening the data and then applying the Euclidean distance ($$\widetilde{x}$$ is the wihtened variable)

  $$d(x,y)=(x-y)^{T}\Sigma^{-1}(x-y)$$


  $$\widetilde{x} = \Sigma^{-1/2}(x-\mu)$$

  $$(\tilde{x_{i}}-\tilde{x_{j}})^{T}(\tilde{x_{i}}-\tilde{x_{j}})=(x_{i}-x_{j})^{T}\Sigma^{-1}(x_{i}-x_{j})$$

- Metric learning is to learn a distance function tuned to a particular task
  - Input: instances of data, and some similarity information
  - Output: a metric function, may in form of

  $$d(x,y)=||f(x)-f(y)||_{2}$$
    
    - The learnt transformation $$f$$ can be linear projection (yield a form of distance same as  Mahalanobis Distance) or non-linear function

  
## Metric Learning with Neural Network

### Classification loss vs. Metric learning loss
- In cases where a classification loss is applicable, why are embeddings used during test time, instead of the logits or the subsequent softmax values? Typically, embeddings are preferred when the task is some variant of information retrieval, where the goal is to return data that is most similar to a query.
- open set classification 
- face verification
- person re-identification
- some time the label is not available, while similarity information is available

### Zeros shot learning 
- <https://github.com/mbsariyildiz/zsl_eval>


### Reference
- <https://slazebni.cs.illinois.edu/spring17/lec09_similarity.pdf>
- 2013, FTML, [Metric Learning: A Survey](https://ieeexplore.ieee.org/document/8186753)
- [A Metric Learning Reality Check](https://arxiv.org/pdf/2003.08505.pdf)
  - <https://github.com/KevinMusgrave/powerful-benchmarker>

 
- triplet network
- siamese network
- deep hash 




- How to generate pair / triplet for training ? 
