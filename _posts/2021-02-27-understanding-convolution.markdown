---
layout: post
title:  "Understanding Convolution"
date:   2021-02-27 12:33:15 +0800
usemathjax: true
categories: jekyll update
---


## Definition

- Continuous case: 

  $$(f*g)(t) := \int_{-\infty}^{\infty}f(\tau)g(t-\tau)d\tau = \int_{-\infty}^{\infty}f(t-\tau)g(\tau)d\tau$$

- Discrete case

   $$(f*g)[n]:=\sum_{m=-\infty}^{\infty}f[m]g[n-m]=\sum_{m=-\infty}^{\infty}f[n-m]g[m]$$

- Relation to polynomial multiplication

  - We have:

    $$
    \begin{align*}
     &(1+2x+3x^2)(3+2x+2x^2) \\
    =& 3(1+2x+3x^2) + 2x(1+2x+3x^2) + 2x^2(1+2x+3x^2) \\
    =& 3 + (3*2+2)x + (3*3+2*2+2)x^2 + (2*3+2*2)x^3 + 2*3x^4 \\
    =& 3 + 8x + 15x^2 + 10x^3 + 6x^4
    \end{align*}
    $$

  - ```python
    np.convolve([1,2,3],[3,2,2]) #array([ 3,  8, 15, 10,  6])
    ```

  - Polynomial multiplication is equivalent to convoluting their coefficients ...

- Relation to  cross correlation

  - Cross correlation is defined as 

    $$(f*g)(t):=\int_{-\infty}^{\infty}\overline{f(\tau)}g(t+\tau)d\tau \\$$

  - $$\overline{f(t)}$$ is complex conjugate of $$f(t)$$
  
  - Cross correlation of $$f(t)$$ and $$g(t)$$ is equivalent to convolution of $$\overline{f(-t)}$$ and $$g(t)$$

## Implementations
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

## Convolution in deep learning
- Convolutional neural network is actually cross correlational neural network. 
- As the convolution kernel is trainable, convolution and cross correlation is equivalent under this context
- See <https://pytorch.org/docs/stable/generated/torch.nn.Conv2d.html>
- There is also a concept called transpose convolution, some one called it "deconvolution" (but it is not really a [deconvolution](https://en.wikipedia.org/wiki/Deconvolution))
- For convolutional network, it's important to understand the shape of input and output feature map. The following material give a clear description
  - <https://github.com/vdumoulin/conv_arithmetic>
  - <https://arxiv.org/abs/1603.07285>
- Convolution and Transpose Convolution is recurrently utilized in architechture like [U-net](https://lmb.informatik.uni-freiburg.de/people/ronneber/u-net/), see this implementation
  - <https://github.com/milesial/Pytorch-UNet/blob/master/unet/unet_parts.py>

## Reference

1. <https://en.wikipedia.org/wiki/Convolution>
2. <https://en.wikipedia.org/wiki/Cross-correlation>
3. <https://numpy.org/doc/stable/reference/generated/numpy.convolve.html>
4. <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.convolve.html>
5. <http://home.cse.ust.hk/~dekai/271/notes/L03/L03.pdf> polynomial multiplication and convolution
6. <https://arxiv.org/pdf/1603.07285.pdf> A guide to convolution arithmetic for deep learning