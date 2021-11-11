---
layout: post
title:  "Manipulate Binary Data"
date:   2021-11-11 20:36:17 +0800
usemathjax: true
categories: jekyll update
---

### In native python


### In numpy

- [packbits](https://numpy.org/doc/stable/reference/generated/numpy.packbits.html): pack binary data into integer (np.uint8) by byte

```python
import numpy as np
X = np.random.randint(low=0,high=2,size=(1000,32)) # generate some binary data
X_packed8 = np.packbits(X,axis=1) # pack every 8 bit by axis 1
# X_packed8.shape
# (1000, 4)
X_packed32 = X_packed8.view(np.uint32)
# X_packed32.shape
# (1000, 1)
```

- [tobyte](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.tobytes.html): save numpy data to byte string

```python
x = np.array([0,1,2,4,8,16,32,64,128,256,512,1024],dtype=np.int16)
x.tobytes()
# b'\x00\x00\x01\x00\x02\x00\x04\x00\x08\x00\x10\x00 \x00@\x00\x80\x00\x00\x01\x00\x02\x00\x04'
```

- [frombuffer](https://numpy.org/doc/stable/reference/generated/numpy.frombuffer.html): load data from byte string