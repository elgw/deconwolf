This folder contains the source documents for the github pages.

The github pages lives in the orphaned branch gh-pages.

To build, commands like these are needed:

``` shell
deactivate # in a virtual environment is in use
python -m venv .venv
source .venv/bin/activate
python -m pip install sphinx sphinx_rtd_theme
make html
```

To build directly to the gh-pages branch one option is to create a
symbolic link to a point where the gh-pages branch is cloned.

``` shell
ln -s ../../../deconwolf-gh-pages/ build/html
```


`gh-pages/source/man` contains the source for the man-pages.

Here are some pointers for how this works:

[https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-man_pages]

[https://github.com/sphinx-doc/sphinx/blob/25d4ae578b759fa7221c826d14dadb6a12994b25/doc/man/sphinx-quickstart.rst?plain=1#L4]
[https://github.com/sphinx-doc/sphinx/blob/25d4ae578b759fa7221c826d14dadb6a12994b25/doc/conf.py]
