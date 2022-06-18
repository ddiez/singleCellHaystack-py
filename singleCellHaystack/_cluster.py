def cluster_genes(adata, haystack_result, method="kmeans", n_clusters=None, random_state=None):
  from sklearn.cluster import KMeans
  import numpy as np
  import pandas as pd
  from ._haystack import calculate_P_dist
  import scipy.stats as ss

  #grid_points = r["info"]["grid.points"]
  #grid_dist = calculate_dist_to_cells(coord, grid_points)
  #density = calculate_density(grid_dist)
  density = haystack_result["info"]["grid_density"]
  ngrid_points = density.shape[1]

  expression = adata.X.toarray()
  ngenes = expression.shape[1]

  # FIXME: vectorize this computation.
  scores = np.zeros([ngenes, ngrid_points])
  for k in range(ngenes):
    scores[k, :] = calculate_P_dist(density, expression[:, k])

  if method=="kmeans":
    res = KMeans(n_clusters=n_clusters, random_state=random_state).fit(scores)
    clusters = res.labels_

  res = pd.DataFrame({
    "gene": adata.var_names,
    "cluster": clusters
  })

  return res

def plot_gene_clusters(adata, gene_clusters, basis="umap", ncols=4, figsize=None, color_map="coolwarm"):
  import numpy as np
  import matplotlib.pyplot as plt

  if basis in adata.obsm.keys():
    basis_key = basis
  elif f"X_{basis}" in adata.obsm.keys():
    basis_key = f"X_{basis}"

  coord = adata.obsm[basis_key]

  clusters = gene_clusters.cluster.unique()
  scores = np.zeros((adata.n_obs, len(clusters)))
  for k in range(len(clusters)):
    genes = gene_clusters[gene_clusters.cluster == clusters[k]].gene.tolist()
    m = adata[:, genes].X
    scores[:,k] = np.squeeze(np.sum(m, axis=1))

  import matplotlib.pyplot as plt
  nclusters = scores.shape[1]
  nrows = int(np.ceil(nclusters/ncols))

  fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)

  for i in range(nrows):
    for j in range(ncols):
      ax[i, j].grid(False)
      ax[i, j].set_axis_off()

  crow=0
  ccol=0
  for k in range(nclusters):
    im = ax[crow, ccol].scatter(coord[:,0], coord[:,1], s=4, c=scores[:,k], cmap=color_map)
    ax[crow, ccol].set_title("Module: " + str(k))
    ax[crow, ccol].set_axis_on()
    ax[crow, ccol].set_xticks([])
    ax[crow, ccol].set_yticks([])
    ax[crow, ccol].set_xlabel("UMAP1")
    ax[crow, ccol].set_ylabel("UMAP2")
    fig.colorbar(im, ax=ax[crow,ccol])
  #  print("crow: " + str(crow) + ", ccol: " + str(ccol))
    ccol = ccol + 1
    if (ccol > ncols - 1):
      crow+=1
      ccol=0

  fig.tight_layout()

  return scores
