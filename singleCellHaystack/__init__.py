""" 

This is a python implementation of `singleCellHaystack <https://github.com/alexisvdb/singleCellHaystack>`_.

"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("singleCellHaystack")
except PackageNotFoundError:
    # package is not installed
    pass

from anndata import AnnData

from ._haystack import *
from ._plot import *
from ._cluster import *
from ._data import *
from ._tools import *
