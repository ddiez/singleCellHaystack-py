import numpy as np
import pandas as pd
import pytest

da = pytest.importorskip("dask.array")

import singleCellHaystack as hs

rng = np.random.default_rng(0)
ncells, ngenes = 200, 50
X = rng.poisson(1.0, size=(ncells, ngenes)).astype(float)
coord = rng.normal(size=(ncells, 2))


def _haystack_kwargs():
  return dict(ngrid_points=10, n_genes_to_randomize=20, n_randomizations=5, random_state=1, verbose=False)


def test_haystack_dask_array():
  Xd = da.from_array(X, chunks=(50, ngenes))
  coordd = da.from_array(coord, chunks=(50, 2))

  res = hs.haystack(Xd, coordd, **_haystack_kwargs())
  assert res
  assert res.result is not None
  sum = res.top_features()
  assert isinstance(sum, pd.core.frame.DataFrame) is True


def test_haystack_dask_matches_numpy():
  Xd = da.from_array(X, chunks=(50, ngenes))

  res_np = hs.haystack(X, coord, **_haystack_kwargs())
  res_dask = hs.haystack(Xd, coord, **_haystack_kwargs())

  a = res_np.result.sort_values("gene").reset_index(drop=True)
  b = res_dask.result.sort_values("gene").reset_index(drop=True)
  assert np.allclose(a["KLD"], b["KLD"])
  assert np.allclose(a["pval"], b["pval"])


def test_haystack_adata_dask_layer():
  ad = pytest.importorskip("anndata")

  adata = ad.AnnData(X=X)
  adata.obsm["X_pca"] = coord
  adata.layers["dask"] = da.from_array(X, chunks=(50, ngenes))

  res = hs.haystack(adata, "pca", layer="dask", **_haystack_kwargs())
  assert res
  assert res.result is not None
