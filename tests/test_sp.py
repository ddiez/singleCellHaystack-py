import scanpy as sc
import singleCellHaystack as hs

adata=sc.datasets.pbmc3k_processed()
adata_raw=adata.raw.to_adata()

coord=adata_raw.obsm["X_umap"]
exprs=adata_raw.X

foo=hs.haystack_sparse(exprs, coord)

print(foo)
