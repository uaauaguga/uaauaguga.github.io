---
layout: post
title:  "Einstein Summation"
date:   2021-05-13 20:49:02 +0800
usemathjax: true
categories: jekyll update
---


- Einstein summation 

## Some background 
- <https://math.stackexchange.com/questions/1192825/why-use-einstein-summation-notation>
- <https://physics.iitm.ac.in/~PH1010/mkj_Lect_03.pdf>
- <http://dslavsk.sites.luc.edu/courses/phys301/classnotes/einsteinsummationnotation.pdf>

## numpy

- See <https://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.einsum.html>

- Transpose / Permutate axis 

```python
a = np.arange(6).reshape((2,3))
#array([[0, 1, 2],
#       [3, 4, 5]])
a.T
#array([[0, 3],
#       [1, 4],
#       [2, 5]])
np.einsum('ij->ji', a)
#array([[0, 3],
#       [1, 4],
#       [2, 5]])
```

- Take diagnal elements

```python
a = np.arange(25).reshape(5,5)
#array([[ 0,  1,  2,  3,  4],
#       [ 5,  6,  7,  8,  9],
#       [10, 11, 12, 13, 14],
#       [15, 16, 17, 18, 19],
#       [20, 21, 22, 23, 24]])
np.diag(a)
np.einsum('ii->i', a)
np.einsum(a,[0,0],[0])
# array([ 0,  6, 12, 18, 24])
```

- Sum by rows

```python
a.sum(axis=1)
np.einsum('ij->i', a)
np.einsum(a,[0,1],[0])
```

- Sum by columns

```python
a.sum(axis=0)
np.einsum('ij->j', a)
np.einsum(a, [1,0], [0])
```

- Matrix production

```python
a = np.arange(6).reshape((2,3))
b = np.arange(12).reshape((3,4))
np.dot(a,b)
np.einsum('ij,jk->ik', a,b)


a = a.T
np.dot(a.T,b)
np.einsum('ji,jk->ik', a,b)
```


## pytorch
- See <https://pytorch.org/docs/stable/generated/torch.einsum.html>
- Also see <https://github.com/arogozhnikov/einops>