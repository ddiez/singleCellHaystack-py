extensions = [
  "sphinx_rtd_theme",
  "sphinx.ext.autodoc",
  "sphinx.ext.autosummary",
  "nbsphinx",
  "myst_parser"
]

#import sphinx_rtd_theme
html_theme = "sphinx_rtd_theme"

project="singleCellHaystack"

html_context = dict(
    display_github=True,      # Integrate GitHub
    github_user='ddiez',   # Username
    github_repo='singleCellHaystack-py',     # Repo name
    github_version='main',  # Version
    conf_py_path='/',    # Path in the checkout to the docs root
)