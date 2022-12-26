import singleCellHaystack as hs

def test_toy():
  adata = hs.load_toy()
  res = hs.haystack(adata, "tsne")
  assert res 
