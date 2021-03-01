---
layout: post
title:  "Manipulate Large Sparse Matrix"
date:   2021-03-01 15:44:40 +0800
usemathjax: false
categories: jekyll update
---

## Manipulate Large Sparse Matrix

### Sparse matrix in scipy
- Mutiple type of sparse matrix is available in scipy
- <https://docs.scipy.org/doc/scipy/reference/sparse.html>
- Quote here
> 1. csc_matrix: Compressed Sparse Column format
> 2. csr_matrix: Compressed Sparse Row format
> 3. bsr_matrix: Block Sparse Row format
> 4. lil_matrix: List of Lists format
> 5. dok_matrix: Dictionary of Keys format
> 6. coo_matrix: COOrdinate format (aka IJV, triplet format)
> 7. dia_matrix: DIAgonal format
> 

- Useful notes

- dok_matrix and lil_matrix can be construct efficiently, COO is also OK
- > it is strongly discouraged to use NumPy functions directly on these matrices because NumPy may not properly convert them for computations, leading to unexpected (and incorrect) results. If you do want to apply a NumPy function to these matrices, first check if SciPy has its own implementation for the given sparse matrix class, or convert the sparse matrix to a NumPy array (e.g., using the toarray() method of the class) first before applying the method.

- > To perform manipulations such as multiplication or inversion, first convert the matrix to either **CSC** or **CSR** format. The lil_matrix format is row-based, so conversion to CSR is efficient, whereas conversion to CSC is less so.

- > All conversions among the CSR, CSC, and COO formats are efficient, linear-time operations.

```python
# A coo matrix example given by scipy documentation
import numpy as np
from scipy.sparse import coo_matrix
row  = np.array([0, 0, 1, 3, 1, 0, 0])
col  = np.array([0, 2, 1, 3, 1, 0, 0])
data = np.array([1, 1, 1, 1, 1, 1, 1])
coo = coo_matrix((data, (row, col)), shape=(4, 4))
# Duplicate indices are maintained until implicitly or explicitly summed
np.max(coo.data)
# 1
coo.toarray()
#array([[3, 0, 1, 0],
#       [0, 2, 0, 0],
#       [0, 0, 0, 0],
#       [0, 0, 0, 1]])
```


### Sparse matrix in pandas
#### Convert a dense dataframe to sparse dataframe

```python
import numpy as np
import pandas as pd

# Generate a dense matrix/dataframe filled with random value
# Set most of the values to nan
df = pd.DataFrame(np.random.randn(10000, 4))
df.iloc[:9998] = np.nan

# Convert dense matrix to sparse matrix
sdf = df.astype(pd.SparseDtype("float", np.nan))

# sdf.dtypes show as follows:
#0    Sparse[float64, nan]
#1    Sparse[float64, nan]
#2    Sparse[float64, nan]
#3    Sparse[float64, nan]
#dtype: object

'dense : {:0.2f} bytes'.format(df.memory_usage().sum() / 1e3) 
# dense : 320.13 bytes
'sparse: {:0.2f} bytes'.format(sdf.memory_usage().sum() / 1e3) 
# sparse: 0.22 bytes
```

#### Convert scipy sparse matrix to sparse dataframe
- Note input to `pd.DataFrame.sparse.from_spmatrix` must be convertable to to `csc_matrix`

```python
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
row  = np.array([0, 3, 1, 0])
col  = np.array([0, 3, 1, 2])
data = np.array([4, 5, 7, 9])
coom = coo_matrix((data, (row, col)), shape=(4, 4))
# Input to pd.DataFrame.sparse.from_spmatrix must have a .tocsc method
# coom.tocsc() is implicitly called
sdf = pd.DataFrame.sparse.from_spmatrix(coom)
```


#### A practical example
- Reference: <https://stackoverflow.com/questions/31661604/efficiently-create-sparse-pivot-tables-in-pandas>

```python
# Import required package
import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype
from scipy.sparse import coo_matrix

# Generate demo data
# Mimic a very sparse matrix, eg, results of single single cell sequencing
X = np.random.randn(100000,10)
X[X<2] = np.nan
index = [str(i+1) for i in range(100000)]
columns = [str(i+1) for i in range(10)]
records = pd.DataFrame(X,index=index,columns=columns).stack().to_frame().to_records()
table = pd.DataFrame.from_records(records)
table.columns = ["gene","sample_id","values"]

# Convert the table into sparse matrix instead of performing pivoting

sample_ids = sorted(table['sample_id'].unique())
genes = sorted(table['gene'].unique())

sample_id_c = CategoricalDtype(sample_ids, ordered=True)
gene_c = CategoricalDtype(genes, ordered=True)

# Assign a integer categrical code to each gene and each sample
# Perform such mapping for every single record
rows = table['gene'].astype(gene_c).cat.codes
columns = table['sample_id'].astype(sample_id_c).cat.codes

# Convert records to coo matrix
sparse_matrix = coo_matrix((table['values'], (rows,columns)),shape=(gene_c.categories.size,sample_id_c.categories.size))

# Convert coo matrix to sparse dataframe
matrix = pd.DataFrame.sparse.from_spmatrix(sparse_matrix,index=genes,columns=sample_ids).astype(pd.SparseDtype("float", np.nan))
```

### Manipulate sparse matrix in R
- Based on [Matrix](https://cran.r-project.org/web/packages/Matrix/Matrix.pdf) package
- A tutorial on Matrix package <https://cran.r-project.org/web/packages/Matrix/vignettes/Intro2Matrix.pdf>
- Matrix is not in R base, the priority is "recommended", same as MASS package
#### Convert a condense Matrix to sparse Matrix

```R
library(Matrix)
M <- Matrix(10 + 1:28, 4, 7)
(M2 <- cBind(-1, M))
M2[, c(2,4:6)] <- 0
M2[2, ] <- 0
M2 <- rBind(0, M2, 0)
M2[1:2,2] <- M2[3,4:5] <- NA
sM <- as(M2, "sparseMatrix")
```

- `dgTMatrix` in `Matrix` is similar to `coo_matrix` in `scipy` 

