import singleCellHaystack as hs
import pandas as pd
import numpy as np

def test_toy():
  adata = hs.load_toy()
  res = hs.haystack(adata, "tsne", random_state=1)
  assert res
  assert res.result is not None
  sum = res.top_features()
  assert isinstance(sum, pd.core.frame.DataFrame) is True
  #sum = res["results"]
  #assert np.all(sum.gene.head(3) == ["gene_62", "gene_79", "gene_339"])
  #assert np.all(sum.KLD.head(3) == [2.092214, 2.308899, 1.840823])

def test_scale_coord():
  adata = hs.load_toy()
  res = hs.haystack(adata, "tsne", scale_coord=False, random_state=1)
  assert res
  assert res.result is not None
  sum = res.top_features()
  assert isinstance(sum, pd.core.frame.DataFrame) is True

def test_kld_new():
  adata = hs.load_toy()
  res = hs.haystack(adata, "tsne", random_state=1, kld_method="new")
  assert res
  assert res.result is not None
  sum = res.top_features()
  assert isinstance(sum, pd.core.frame.DataFrame) is True
  #sum = res["results"]
  #assert np.all(sum.gene.head(3) == ["gene_62", "gene_79", "gene_339"])
  #assert np.all(sum.KLD.head(3) == [2.092214, 2.308899, 1.840823])

def test_HaystackResult():
  r = pd.DataFrame({"gene": ["A", "B"], "logpval_adj": [0.7, 0.01]})
  res = hs.HaystackResult(result=r, info="goo")
  sum = res.top_features()
  #assert sum
  assert isinstance(sum, pd.core.frame.DataFrame) is True