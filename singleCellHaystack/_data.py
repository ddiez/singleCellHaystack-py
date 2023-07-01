import importlib.resources

def load_toy():
    import scanpy as sc
    path = importlib.resources.files(__name__).joinpath("datasets/toy.h5ad")
    return sc.read_h5ad(path)
