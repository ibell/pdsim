# Jupyter notebooks

# Setup

Jupyter notebooks allow for interactive coding in Python, R, Julia, and many other programming languages in a cell-based architecture. The easiest way to get access to the Jupyter notebook is to install it via conda. If you install Anaconda, jupyter and its dependent packages are included directly.  Otherwise, you can install jupyter notebook from conda:

```
conda install jupyter notebook
```

No matter how you have obtained jupyter notebook, you can open a shell in the folder where you would like to start jupyter notebook and do

```
jupyter notebook
```
and your browser should start and display the home screen.  You can click on folders to go into them, or on files with an ``.ipynb`` file extension to open them for editing and execution

# Notes

* Shift+Enter executes a cell and moves to the next cell
* Ctrl+Enter executes a cell and stays in the same cell
* You can also change the cell type to Markdown and add Markdown to it (which is a superset of HTML, so HTML is also allowed).  LaTeX is also supported in Markdown.  Some information from github: https://guides.github.com/features/mastering-markdown/
* Some more information here: https://realpython.com/jupyter-notebook-introduction/