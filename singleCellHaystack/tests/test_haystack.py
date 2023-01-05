import singleCellHaystack as hs
import pandas as pd

def test_toy():
  adata = hs.load_toy()
  res = hs.haystack(adata, "tsne")
  assert res
  assert res["results"] is not None
  assert isinstance(res["results"], pd.core.frame.DataFrame) is True
