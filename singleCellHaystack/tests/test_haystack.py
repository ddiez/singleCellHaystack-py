import singleCellHaystack as hs
import pandas as pd
import numpy as np

def test_toy():
  adata = hs.load_toy()
  res = hs.haystack(adata, "tsne", random_state=1)
  assert res
  assert res["results"] is not None
  assert isinstance(res["results"], pd.core.frame.DataFrame) is True
  #sum = res["results"]
  #assert np.all(sum.gene.head(3) == ["gene_62", "gene_79", "gene_339"])
  #assert np.all(sum.KLD.head(3) == [2.092214, 2.308899, 1.840823])

def test_kld_new():
  adata = hs.load_toy()
  res = hs.haystack(adata, "tsne", random_state=1, kld_method="new")
  assert res
  assert res["results"] is not None
  assert isinstance(res["results"], pd.core.frame.DataFrame) is True
  #sum = res["results"]
  #assert np.all(sum.gene.head(3) == ["gene_62", "gene_79", "gene_339"])
  #assert np.all(sum.KLD.head(3) == [2.092214, 2.308899, 1.840823])
