#!/bin/sh

CONDA_DIR=~/miniconda3
VERSION=0.1.1.dev0
PACKAGE=franklab_mountainsort

conda build conda.recipe/meta.yaml -c conda-forge -c flatiron -c edeno  --no-anaconda-upload
conda convert $CONDA_DIR/conda-bld/osx-64/$PACKAGE-$VERSION*.tar.bz2 -o $CONDA_DIR/conda-bld/ --platform all
anaconda upload $CONDA_DIR/conda-bld/*/$PACKAGE-$VERSION*.tar.bz2
