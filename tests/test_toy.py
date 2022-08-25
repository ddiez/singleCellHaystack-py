import singleCellHaystack as hs
adata = hs.load_toy()
res = hs.haystack(adata, coord="tsne")
print(res["results"])
