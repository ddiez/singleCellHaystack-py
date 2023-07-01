import importlib_resources

def load_toy():
    import scanpy as sc
    path = importlib_resources.files(__name__).joinpath("datasets/toy.h5ad")
    return sc.read_h5ad(path)
