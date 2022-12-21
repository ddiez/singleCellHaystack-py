singleCellHaystack-py
=====================

[![Lifecycle:beta](https://img.shields.io/badge/lifecycle-beta-orange.svg)]()
[![PyPI](https://img.shields.io/pypi/v/singleCellHaystack?logo=PyPI)](https://pypi.org/project/singleCellHaystack)
[![PyPIDownloads](https://pepy.tech/badge/singleCellHaystack)](https://pepy.tech/project/singleCellHaystack)

This repository contains a python implementation of [singleCellHaystack](https://github.com/alexisvdb/singleCellHaystack).

# Installation

For now, clone this repository, then run the following from the repository folder:

```
pip install singleCellHaystack
```

# Example

```{python}
import scanpy as sc
import singleCellHaystack as hs

adata = sc.read_h5ad("data.h5ad")

res = hs.haystack(adata)
```
