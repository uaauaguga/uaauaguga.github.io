---
layout: post
title:  "Convolutions In Python"
date:   2021-02-27 12:33:15 +0800
usemathjax: true
categories: jekyll update
---

## Convolution Operation in Python

### Definition

- Continuous case: 
  - $$f(t*g)(t) := \int_{-\infty}^{\infty}f(\tau)g(t-\tau)d\tau = \int_{-\infty}^{\infty}f(t-\tau)g(\tau)d\tau$$
- Discrete case
  - $$(f*g)[n]=\sum_{m=-\infty}^{\infty}f[m]g[n-m]$$
- Relation to polynomial multiplication
- Relation to  cross correlation

### Implementations
- `numpy.convolve`

  ```python
  ## https://numpy.org/doc/stable/reference/generated/numpy.convolve.html
  np.convolve([1, 2, 3], [0, 1, 0.5],"full") #array([0. , 1. , 2.5, 4. , 1.5]) N+M-1
  np.convolve([1, 2, 3], [0, 1, 0.5],"same") #array([1. , 2.5, 4. ]) max(M, N)
  np.convolve([1, 2, 3], [0, 1, 0.5],"valid") #array([2.5]) max(M, N) - min(M, N) + 1
  ```

- `scipy.signal.convolve`

  - Similar to numpy

- `scipy.signal.convolve2d`

- `scipy.signal.fftconvolve`

- A naive toy implementation for the "valid" mode convolution and cross correlation

  - Reverse the kernel, sliding and get inner products 

  ```python
  def convolution(x,y):
      z = []
      x,y = np.array(x),np.array(y)
      y = y[::-1]
      for i in range(0,x.shape[0] - y.shape[0]+1):
          z.append(np.dot(x[i:i+y.shape[0]],y))
      return np.array(z)
          
  def crosscorrelation(x,y):
      z = []
      x,y = np.array(x),np.array(y)
      for i in range(0,x.shape[0] - y.shape[0]+1):
          z.append(np.dot(x[i:i+y.shape[0]],y))
      return np.array(z)
  ```

### Applications

- Smoothing with a moving average kernel
- Image processing