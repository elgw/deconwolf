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
