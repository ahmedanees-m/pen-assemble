"""Sphinx configuration for PEN-ASSEMBLE documentation."""
import sys
from pathlib import Path

# Add repository root to path so autodoc can import pen_assemble
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# -- Project information -------------------------------------------------------
project   = "PEN-ASSEMBLE"
copyright = "2026, Anees Ahmed"
author    = "Anees Ahmed"
release   = "0.5.0"
version   = "0.5"

# -- General configuration -----------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
]

templates_path   = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Napoleon settings (Google / NumPy docstrings)
napoleon_google_docstring  = True
napoleon_numpy_docstring   = True
napoleon_include_init_with_doc = True

# autodoc
autodoc_member_order        = "bysource"
autodoc_typehints           = "description"
autodoc_typehints_format    = "short"
always_document_param_types = True

# intersphinx
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "numpy":  ("https://numpy.org/doc/stable", None),
}

# -- HTML output ---------------------------------------------------------------
html_theme = "furo"
html_title = "PEN-ASSEMBLE v0.5.0"
html_static_path = ["_static"]
html_theme_options = {
    "light_logo": "pen_assemble_logo.png",
    "dark_logo":  "pen_assemble_logo.png",
}

# -- doctest -------------------------------------------------------------------
doctest_global_setup = "import pen_assemble"
