extensions = [
  "sphinx_book_theme",
  "sphinx.ext.autodoc",
  "sphinx.ext.autosummary",
  "nbsphinx",
  "myst_parser"
]

#import sphinx_rtd_theme
html_theme = "sphinx_book_theme"

project="singleCellHaystack"

html_theme_options = {
  "path_to_docs": "docs/",
  "use_repository_button": True,
  "repository_url": "https://github.com/ddiez/singleCellHaystack-py"
}