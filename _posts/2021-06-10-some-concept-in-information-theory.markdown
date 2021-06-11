---
layout: post
title:  "Basic Concepts in Information Theory"
date:   2021-06-10 14:14:50 +0800
usemathjax: true
categories: jekyll update
---

### Information content
- 越稀少的事件带来的surprise越大

$$I(X) = -ln(p_X(x)) = ln \frac{1}{p_X(x)}$$

### Entropy
- 熵是信息量的期望

- 离散的情形

$$H(X) = \sum_{i}p_{i}ln(p_{i})$$

- 连续的情形下叫做微分熵

$$H(X) \equiv \mathbb{E}_{X}I(X) = -\int_{-\infty}^{\infty}p_X(x)ln(p_X(x))dx$$

- a gaussian variable has the largest entropy among all random variables of equal variance

- 条件分布的熵

- The conditional entropy $$H(Y \mid X)$$ is the average additional informaion needed to specified `Y`

$$p_{X,Y}(x,y) = p_{Y|X}(y|x)p_{X}(x)$$

$$p_{X,Y}(x,y)ln(p_{X,Y}(x,y)) = p_{X,Y}(x,y)ln(p_{Y|X}(y|x)) + p_{X,Y}(x,y)ln(p_{X}(x))$$


$$H(Y|X) = \int_{-\infty}^{\infty}p_{X}(x)H(Y|X=x)dx = \int_{-\infty}^{\infty} p_{X}(x) \int_{-\infty}^{\infty} p_{Y|X=x}(y)ln(p_{Y|X=x}(y))dydx$$

$$\int_{-\infty}^{\infty}\int_{-\infty}^{\infty} p_{X,Y}(x,y)ln(p_{X,Y}(x,y))dxdy = \int_{-\infty}^{\infty}\int_{-\infty}^{\infty} p_{Y|X=x}(y)p_{X}(x)ln(p_{Y|X=x}(y))dxdy + \int_{-\infty}^{\infty}\int_{-\infty}^{\infty} p_{Y|X=x}(y)p_{X}(x)ln(p_{X}(x))dxdy$$


$$H((X,Y)) = H(Y|X) + H(X)$$

- The bayes rule for conditional entropy

$$H((X,Y)) \equiv H(X) + H(Y|X) \equiv H(X|Y) + H(Y)$$

$$H(X|Y)  = H(X) + H(Y|X) - H(Y)$$


### Cross entropy
- 交叉熵

$$H(X,Y) \equiv \mathbb{E}_{X}I(Y) = -\int_{-\infty}^{\infty}p(x)ln(q(x))dx$$

- 交叉熵和分类问题中的负对数似然损失函数是等价的，见<https://en.wikipedia.org/wiki/Cross_entropy>

- 注意交叉熵和联合分布的熵不是一个东西！<https://math.stackexchange.com/questions/2505015/relation-between-cross-entropy-and-joint-entropy>

### KL divergence
- KL散度/相对熵

$$D_{KL}(X\|Y) \equiv H(X,Y) - H(Y) = \int_{-\infty}^{\infty}p(x)ln \frac{q(x)}{p(x)}dx= -\int_{-\infty}^{\infty}p(x)ln \frac{p(x)}{q(x)}dx$$

- Note $$D_{KL}(X\|Y) \neq D_{KL}(Y\|X)$$


### Mutual information
- 互信息/互传信息量
- 如果$$X$$和$$Y$$独立，我们有 $$p_{x,y}(X,Y) = p_{X}(x)p_{Y}(y)$$
- 互信息$$X$$和$$Y$$的互信息定义为联合分布 $$(X,Y)$$ 和 服从$$p_{X}(x)p_{Y}(y)$$的随机变量之间的KL divergence 

$$MI(X,Y) \equiv D_{KL}(p_{X,Y}\|p_{X}p_{Y}) = \int_{-\infty}^{\infty} p_{X,Y}(x,y)ln \frac{p_{X,Y}(x,y)}{p_{X}(x)p_{Y}(y)}$$

- Some relationship 

$$MI(X,Y) \equiv H(X) - H(X|Y) \equiv H(Y) - H(Y|X) \equiv H((X,Y)) - H(X|Y) - H(Y|X)$$

### Reference

- <https://www2.isye.gatech.edu/~yxie77/ece587/Lecture2.pdf>
- <https://en.wikipedia.org/wiki/Kullback–Leibler_divergence>
