import scanpy as sc
import singleCellHaystack as hs

adata=sc.datasets.pbmc3k_processed()
adata_raw=adata.raw.to_adata()

coord=adata_raw.obsm["X_umap"]
exprs=adata_raw.X
genes=adata_raw.var_names

foo=hs.haystack_array(exprs, coord, features=genes)
print(foo["results"])

foo2=hs.haystack_adata(adata_raw, "X_umap")
print(foo2["results"])

foo3=hs.haystack(adata_raw, "X_umap")
print(foo3["results"])

foo4=hs.haystack(exprs, coord, features=genes)
print(foo4["results"])

import numpy as np
foo5=hs.haystack(np.asarray(exprs.todense()), coord, features=genes)
print(foo5["results"])
