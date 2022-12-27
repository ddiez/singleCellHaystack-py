singleCellHaystack: finding neddles (genes) in haystacks (single cell genomics data)
====================================================================================

[![Lifecycle:beta](https://img.shields.io/badge/lifecycle-beta-orange.svg)](https://github.com/ddiez/singleCellHaystack-py)
[![](https://github.com/ddiez/singleCellHaystack-py/actions/workflows/python-package.yml/badge.svg)](https://github.com/ddiez/singleCellHaystack-py/actions/workflows/python-package.yml)
[![PyPI](https://img.shields.io/pypi/v/singleCellHaystack?logo=PyPI)](https://pypi.org/project/singleCellHaystack)
[![PyPIDownloads](https://pepy.tech/badge/singleCellHaystack)](https://pepy.tech/project/singleCellHaystack)

This repository contains a python implementation of [singleCellHaystack](https://github.com/alexisvdb/singleCellHaystack).

# Installation

You can install singleCellHaystack from [pypi](https://pypi.org):

```
pip install singleCellHaystack
```

# Example

```
import scanpy as sc
import singleCellHaystack as hs

adata = sc.read_h5ad("data.h5ad")

res = hs.haystack(adata)
```
