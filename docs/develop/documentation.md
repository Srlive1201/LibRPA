# Documentation

## Overview

LibRPA uses [Sphinx](https://www.sphinx-doc.org/en/master/) toolchain to build its user documentation.
C++ API documentation is generated using [Doxygen](https://www.doxygen.nl)
and connected to Sphinx using [Doxysphinx](https://github.com/boschglobal/doxysphinx)
and [doxylink](https://sphinxcontrib-doxylink.readthedocs.io/en/stable/).
We also use the following packages to make the writing easier and prettify the documentation
<!-- - [Sphinx book theme](https://github.com/executablebooks/sphinx-book-theme) -->
- [Sphinx RTD theme](https://github.com/readthedocs/sphinx_rtd_theme)
- [Sphinx design](https://github.com/executablebooks/sphinx-design)
- [Sphinx copy-button](https://github.com/executablebooks/sphinx-copybutton)
- [MyST parser](https://github.com/executablebooks/MyST-Parser)
- [Doxygen Awesome](https://jothepro.github.io/doxygen-awesome-css/)

## Build

To build the documentation, the mentioned prerequisites can be installed from PyPi
```shell
pip3 install -r docs/requirements.txt
pip3 install doxysphinx
```
or using conda
```shell
conda install -c conda-forge --file docs/requirements.txt
# unfortunately doxysphinx can only be installed from pypi
pip3 install doxysphinx
```
It is recommened to create a dedivated virtual environment to install the packages.

After successful install of the prerequisites, the document can be built by issuing
the following `make` command under `docs` directory.
```shell
make html
# or simply
make
```
This will download the required style sheet and run the toolchain.
The generated (HTML) documentation is under `_build/html` directory.
You can preview the documentation by opening `_build/html/index.html` using any
web browser.
