singleCellHaystack
==================

[![Lifecycle:beta](https://img.shields.io/badge/lifecycle-beta-orange.svg)](https://github.com/ddiez/singleCellHaystack-py)
[![](https://github.com/ddiez/singleCellHaystack-py/actions/workflows/python-package.yml/badge.svg)](https://github.com/ddiez/singleCellHaystack-py/actions/workflows/python-package.yml)
[![PyPI](https://img.shields.io/pypi/v/singleCellHaystack?logo=PyPI)](https://pypi.org/project/singleCellHaystack)
[![PyPIDownloads](https://pepy.tech/badge/singleCellHaystack)](https://pepy.tech/project/singleCellHaystack)

This repository contains a python implementation of [singleCellHaystack](https://github.com/alexisvdb/singleCellHaystack) (version >= 1.0.0).

This package is currently in beta. The most important functionality in the R package works, but some features are not yet available. Here is a (probably imcomplete) list of missing features. Some will be added in the future.

* `weights.advanced.Q` (formerly known as `use.advanced.sampling`).
* `seeding` method for calculating grid points.

# Installation

You can install singleCellHaystack from [pypi](https://pypi.org):

```
pip install singleCellHaystack
```

>Support for conda installation will be added in the future.

# Example

```
import scanpy as sc
import singleCellHaystack as hs

adata = sc.read_h5ad("data.h5ad")

[... process adata object ...]

res = hs.haystack(adata, coord="pca")
res["results"]
```

# References

- Our manuscript describing the updated, more generally applicable version of `singleCellHaystack` inclusing this python implementation is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.11.13.516355v1).

- Our manuscript describing the original implementation of `singleCellHaystack` for R ([version 0.3.4](https://github.com/alexisvdb/singleCellHaystack/tree/binary)) was published in [Nature Communications](https://doi.org/10.1038/s41467-020-17900-3).

If you use `singleCellHaystack` in your research please cite our work using:

Vandenbon A, Diez D (2020). “A clustering-independent method for finding differentially expressed genes in single-cell transcriptome data.” *Nature Communications*, *11*(1), 4318. [doi:10.1038/s41467-020-17900-3](https://doi.org/10.1038/s41467-020-17900-3).
