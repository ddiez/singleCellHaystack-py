singleCellHaystack-py
=====================

[![Lifecycle:beta](https://img.shields.io/badge/lifecycle-beta-orange.svg)]()

This repository contains a python implementation of [singleCellHaystack](https://github.com/alexisvdb/singleCellHaystack).

# Installation

For now, clone this repository, then run the following from the repository folder:

```
pip install .
```

# Example

```{python}
import scanpy as sc
import singleCellHaystack as hs

adata = sc.read_h5ad("data.h5ad")

res = hs.haystack(adata)
```
