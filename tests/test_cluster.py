import singleCellHaystack as hs
adata = hs.load_toy()
res = hs.haystack(adata, coord="tsne")
gene_mods = hs.cluster_genes(adata, res.head(100), n_clusters=3)
print(gene_mods)
