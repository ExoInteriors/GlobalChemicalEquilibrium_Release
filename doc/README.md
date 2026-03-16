# Documentation

Author of the documentation: Simon Grimm

The documentation is built with sphinx


## How to install sphinx

It is recommended to create a virtual environment to install all necessary packages

### Create a virtual environment

python3 -m venv sphinx-env


### Load the virtual environment

source sphinx-env/bin/activate

### Install the following packages in the virtual environment

sudo apt-get install python3-sphinx

sudo apt-get install latexmk

pip3 install sphinx-astropy

pip3 install astropy

pip3 install sphinxcontrib-bibtex

pip3 install sphinx_rtd_theme


## How to build the documentation

Load the virtual environment.

cd doc


### create a html documentation

make html

the documentation is then build at doc/build/html/index.html

### create a pdf documentation

make latexpdf

the documentation is then build at doc/build/latex/fastglobalchemcode.pdf 

