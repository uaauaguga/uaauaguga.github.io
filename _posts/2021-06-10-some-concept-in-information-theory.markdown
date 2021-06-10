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
- 这种连续的情形有时又叫做微分熵

$$H(X) = \mathbb{E}_{X}I(X) = -\int_{-\infty}^{\infty}p_X(x)ln(p_X(x))dx$$

- 条件分布的熵

$$p_{X,Y}(x,y) = p_{Y|X}(y|x)p_{X}(x)$$

$$p_{X,Y}(x,y)ln(p_{X,Y}(x,y)) = p_{X,Y}(x,y)ln(p_{Y|X}(y|x)) + p_{X,Y}(x,y)ln(p_{X}(x))$$

$$H(Y|X) = \int_{-\infty}^{\infty}p_{X}(x)H(Y|X=x)dx = \int_{-\infty}^{\infty} p_{X}(x) \int_{-\infty}^{\infty} p_{Y|X=x}(y)ln(p_{Y|X=x}(y))dydx$$

$$\int_{-\infty}^{\infty}\int_{-\infty}^{\infty} p_{X,Y}(x,y)ln(p_{X,Y}(x,y))dxdy = \int_{-\infty}^{\infty}\int_{-\infty}^{\infty} p_{Y|X=x}(y)p_{X}(x)ln(p_{Y|X=x}(y))dxdy + \int_{-\infty}^{\infty}\int_{-\infty}^{\infty} p_{Y|X=x}(y)p_{X}(x)ln(p_{X}(x))dxdy$$


$$H((X,Y)) = H(Y|X) + H(X)$$


### Cross entropy
- 交叉熵

$$H(X,Y) = \mathbb{E}_{X}I(Y) = -\int_{-\infty}^{\infty}p(x)ln(q(x))dx$$

- <https://en.wikipedia.org/wiki/Cross_entropy>

### KL divergence
- KL散度/相对熵

$$D_{KL}(X||Y) = H(Y|X) = H((X,Y)) - H(Y) = \int_{-\infty}^{\infty}p(x)ln \frac{q(x)}{p(x)}dx= -\int_{-\infty}^{\infty}p(x)ln \frac{p(x)}{q(x)}dx$$


#### Mutual information
- 互信息/互传信息量

$$I(X;Y) = D_{KL}(p_{X,Y}||p_{X}p_{Y}) = \int_{-\infty}^{\infty} p_{X,Y}(x,y)ln \frac{p_{X,Y}(x,y)}{p_{X}(x)p_{Y}(y)}$$


### Reference

- <https://www2.isye.gatech.edu/~yxie77/ece587/Lecture2.pdf>
