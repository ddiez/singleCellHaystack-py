singleCellHaystack
==================

[![Lifecycle:beta](https://img.shields.io/badge/lifecycle-beta-orange.svg)](https://github.com/ddiez/singleCellHaystack-py)
[![](https://github.com/ddiez/singleCellHaystack-py/actions/workflows/python-package.yml/badge.svg)](https://github.com/ddiez/singleCellHaystack-py/actions/workflows/python-package.yml)
[![PyPI](https://img.shields.io/pypi/v/singleCellHaystack?logo=PyPI)](https://pypi.org/project/singleCellHaystack)
[![Downloads](https://static.pepy.tech/badge/singleCellHaystack)](https://pepy.tech/project/singleCellHaystack)

This repository contains a Python implementation of [singleCellHaystack](https://github.com/alexisvdb/singleCellHaystack) (version >= 1.0.0).

This package is currently in beta. The most important functionality in the R package works, but some features are not yet available. Here is a (probably imcomplete) list of missing features. Some will be added in the future.

* `weights.advanced.Q` (formerly known as `use.advanced.sampling`).
* `seeding` method for calculating grid points.
* Hierarchical clustering method for `cluster_genes`.

# Installation

You can install singleCellHaystack from [PyPI](https://pypi.org):

```
pip install singleCellHaystack
```

You can install singleCellHaystack from GitHub with:

```
pip install git+http://github.com/ddiez/singleCellHaystack-py
```

>Support for conda installation will be added in the future.

# Example

```
import scanpy as sc
import singleCellHaystack as hs

adata = sc.read_h5ad("data.h5ad")

[... process adata object ...]

res = hs.haystack(adata, basis="pca")
res.top_features(n=10)
```

# References

- Our manuscript describing the updated, more generally applicable version of `singleCellHaystack` including this Python implementation was published in [Scientific Reports](https://doi.org/10.1038/s41598-023-38965-2).

- Our manuscript describing the original implementation of `singleCellHaystack` for R ([version 0.3.4](https://github.com/alexisvdb/singleCellHaystack/tree/binary)) was published in [Nature Communications](https://doi.org/10.1038/s41467-020-17900-3).

If you use `singleCellHaystack` in your research please cite our work using:

<p>Vandenbon A, Diez D (2023).
&ldquo;A universal tool for predicting differentially active features in single-cell and spatial genomics data.&rdquo;
<em>Scientific Reports</em>, <b>13</b>(1), 11830.
<a href="https://doi.org/10.1038/s41598-023-38965-2">doi:10.1038/s41598-023-38965-2</a>. 
</p>
