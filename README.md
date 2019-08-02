# franklab_mountainsort

```bash
# Make sure base conda environment is up to date.
conda update -n base -c defaults conda
# Create conda environment `mountainlab` and activate
conda create -n mountainlab
conda activate mountainlab

# Install mountain lab packages and processors
conda install -c flatiron -c conda-forge \
			mountainlab \
			mountainlab_pytools \
			ml_ephys \
			ml_ms3 \
			ml_ms4alg \
			ml_pyms \
			ephys-viz \
			qt-mountainview \
			spikeforestwidgets

# Check to see if processors are installed. Should see more than 1.
ml-list-processors

# Install Frank lab custom processors (future note: make conda installable)
git clone https://bitbucket.org/franklab/franklab_msdrift.git
git clone https://bitbucket.org/franklab/franklab_mstaggedcuration.git

# Symlink files
cd $CONDA_PREFIX/etc/mountainlab/packages
ln -s franklab_msdrift .
ln -s franklab_mstaggedcuration .

# Check symlinks
ls -l
```
