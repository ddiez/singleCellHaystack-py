import singleCellHaystack as hs
adata = hs.load_toy()
res = hs.haystack(adata, basis="tsne")
print(res["results"])
