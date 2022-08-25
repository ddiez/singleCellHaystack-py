import pkg_resources

def load_toy():
    import scanpy as sc
    stream = pkg_resources.resource_stream(__name__, 'datasets/toy.h5ad')
    return sc.read_h5ad(stream)
