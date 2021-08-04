---
layout: post
title:  "Deep Hashing for Similar Item Retrieval"
date:   2021-04-27 23:32:02 +0800
usemathjax: true
categories: jekyll update
---


## Paper list
- <https://github.com/caoyue10/DeepHash-Papers>

## Implementations
- <https://github.com/thulab/DeepHash>
- <https://github.com/thuml/HashNet/tree/master/pytorch#datasets>
- <https://github.com/swuxyj/DeepHash-pytorch>
- <https://github.com/flyingpot/pytorch_deephash>

## Unsupervised hashing techniques
- iterative quantization (ITQ)
- spectral hashing
- discrete graph hashing
- spherical Hashing 
- anchor graph hashing 
- stochastic generative hashing 

## Reading Notes

### Review on deep hash
- 2020, [A Survey on Deep Hashing Methods](https://arxiv.org/abs/2003.03369)

- Network Architecture
  - Shallower architectures (AlexNet, CNN-F) for simple datasets (MINST, CIFAR-10)
  - Deep architectures (VGGNet, ResNet50) for complex datasets such as NUSWIDE and COCO

- Similarity measurement
  - Input space: the ground truth
  - Embedding space and hash space 

- The loss fucntion
  - The similarity loss: preserve the similarity ordering in input space to hash space
  - Classification loss: sematic label of input samples may be utilized
  - Quantization loss: push embedding close to -1 or 1
    - pairwise quantization loss
    - Cauchy quantization loss
  - Other regularizing terms

- Optimization
  - the sign function is non-differentiable, and can not produce back propogate gradient
  - so-called ill-posed gradient problem
  - Several solution
    - add quantization loss as a penalty term
    - only use back propagation for solving a subproblem
    - adopt that continuous relaxation by replacing sign function with **tanh** or **sigmoid**


### Decorrelate the embedding

- We expect compact hash code, but sometimes different bits in the learnt embedding can be highly correlated, lead to information redundancy
- Some study decorrelate the embedding by an additional loss term
- Denote the embedded binary code as $$B$$, then the decorrelation loss can be defined as 

$$L_{DeCov} = \frac{1}{2}(\|C\|_{F}^{2}-\|diag(C)\|_{F}^{2})$$

- $$C$$ is the covariance matrix of the binary embedding  

$$C_{i,j}=\frac{1}{N}\sum_{n}(B_{i}^{n}-\mu_{i})(B_{j}^{n}-\mu_{j})^T$$

- Reference 
  - 2015, CVPR, [Deep Hashing for Compact Binary Codes Learning](https://www.cv-foundation.org/openaccess/content_cvpr_2015/papers/Liong_Deep_Hashing_for_2015_CVPR_paper.pdf)
  - 2015, IJCAI, [Deep Multimodal Hashing with Orthogonal Regularization](https://www.ijcai.org/Proceedings/15/Papers/324.pdf)
  - 2016, ICLR, [Reducing Overfitting in Deep Networks by Decorrelating Representations](https://arxiv.org/pdf/1511.06068.pdf)
  - 2017, AAAI, [Pairwise Relationship Guided Deep Hashing for Cross-Modal Retrieval](https://aaai.org/ocs/index.php/AAAI/AAAI17/paper/view/14326)
  - 2017, ICLR, [Regularizing CNNs with Locally Constrained Decorrelations](https://arxiv.org/abs/1611.01967)
  - 2018, AAAI, [On Trivial Solution and High Correlation Problems in Deep Supervised Hashing](https://eprints.lancs.ac.uk/id/eprint/123575/1/2018_3.pdf)
  - 2018, NIPS, [Can We Gain More from OrthogonalityRegularizations in Training Deep CNNs?](https://arxiv.org/pdf/1810.09102.pdf)

- Also have a look at <https://paperswithcode.com/method/orthogonal-regularization>

### Balancing the code usage


### HashNet
- Subject to the ill-posed gradient difficulty in the optimization with sign activations, existing deep learning to hash methods need to first learn continuous representations and then generate binary hash codes in a separated binarization step, which suffer from substantial loss of retrieval quality


### DPSH
- <https://github.com/TreezzZ/DPSH_PyTorch/blob/master/dpsh.py>
- <https://github.com/jiangqy/DPSH-pytorch>
- <https://github.com/swuxyj/DeepHash-pytorch>

- 2017, NIPS, [Deep Supervised Discrete Hashing](https://arxiv.org/abs/1705.10999)

### Notes

- For binary code $$b_{i}$$ and $$b_{j}$$ of $$K$$ bits, their hamming distance is:

$$dist_{H} = \frac{1}{2} (K-b_{i}^{T}b_{j})$$

- Suppose we have $$N$$ samples, the label contains $$C$$ categories, a binary encoded label is $$Y_{c,N}$$
- According to the label matrix, we have a similarity matrix $$S_{N,N}$$, $$S_{i,j}=1$$ if sample $$i$$ and sample $$j$$ shares same label, otherwise $$S_{i,j}=0$$
- Suppose the binary codes is $$B_{K,N}$$. Given the similarity matrix $$S_{N,N}$$, the posterior probability of the binary code is:


$$p(B|S) \propto P(S|B)P(B) = \sum_{s_{i,j} \in S}{p(s_{i,j}|B)p(B)} $$

- We can see the larger the dot product, the smaller the hamming distance

- Use a sigmoid function for modeling the relation between $$s_{i,j}$$ and learnt binary code

$$p(s_{i,j}=1|B) = \frac{1}{1+e^{-\frac{1}{2}b_{i}^{T}b_{j}}}$$

- We use negative conditional loglikelihood of similarity matrix as dissimilarity loss 

$$
\begin{align*}
   J =& -\ln p(S|B) = - \sum_{s_{i,j} \in S}{p(s_{i,j}|B)}\\
     =& -\sum_{s_{i,j}=0}{p(s_{i,j}|B)}-\sum_{s_{i,j}=1}{p(s_{i,j}|B)} \\
     =& -\sum_{s_{i,j}=1}{\frac{1}{2}b_{i}^{T}b_{j}-\ln(1+e^{\frac{1}{2}b_{i}^{T}b_{j}})}
        -\sum_{s_{i,j}=0}{\ln{1}-\ln(1+e^{\frac{1}{2}b_{i}^{T}b_{j}})} \\
     =& -\sum_{s_{i,j}=1}{\frac{1}{2}b_{i}^{T}b_{j}*1-\ln(1+e^{\frac{1}{2}b_{i}^{T}b_{j}})}
        -\sum_{s_{i,j}=0}{\frac{1}{2}b_{i}^{T}b_{j}*0-\ln(1+e^{\frac{1}{2}b_{i}^{T}b_{j}})} \\
     =& -\sum_{s_{i,j} \in S}{(\frac{1}{2}s_{i,j}b_{i}^{T}b_{j}-\ln(1+e^{\frac{1}{2}b_{i}^{T}b_{j}}))}
\end{align*}
$$

- Ideal binary code should achieve good classification performance

- Suppose $$W_{K,C}$$ is the weight parameter, $$Y=W^{T}B$$

- The classification loss can be written as 


$$Q = \sum_{i=1}^{N}L(y_{i},W^{T}B) + \lambda||W||^2_F$$

- it is essential to keep the discrete nature of the binary codes
- Solution: use a  auxiliary variable

### Implementations

- <https://github.com/TreezzZ/DSDH_PyTorch>
- <https://github.com/liqi-casia/DSDH-HashingCode>

