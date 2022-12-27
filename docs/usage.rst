Usage
=====

To use singleCellHaystack you need to do::

  import singleCellHaystack as hs
  
  adata = hs.load_toy()
  hs.haystack(adata, "tsne")
