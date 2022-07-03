singleCellHaystack-py
=====================

This repository contains a python implementation of [singleCellHaystack](https://github.com/alexisvdb/singleCellHaystack).


```{python}
import scanpy as sc
import singleCellHaystack as hs

adata = sc.read_h5ad("data.h5ad")

res = hs.haystack(adata)
```
