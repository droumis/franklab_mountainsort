# franklab_mountainsort

## Steps to Spike Sort
1. Install
```bash
conda install -c edeno franklab_mountainsort
```
2. Move data using `franklab_mountainsort.move_mda_data`
3. Run spike sorting `franklab_mountainsort.run_spike_sorting`

## Old
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

# Check processors
ml-list-processors

```

You should see at least the following processors:

**ml_ephys**
+ ephys.bandpass_filter
+ ephys.synthesize_random_waveforms
+ ephys.synthesize_random_firings
+ ephys.synthesize_timeseries
+ ephys.whiten

**ml_ms3**
+ ms3.cluster_metrics
+ ms3.combine_cluster_metrics
+ ms3.isolation_metrics
+ ms3.mask_out_artifacts

**ml_ms4**
+ ms4alg.apply_label_map
+ ms4alg.create_label_map
+ ms4alg.sort

**ml_pyms**
+ pyms.extract_timeseries
+ pyms.extract_clips

**franklab_mstaggedcuration**
+ pyms.add_curation_tags
+ pyms.merge_burst_parents

**franklab_msdrift**
+ pyms.anneal_segments


### Next Move data

ml-link-python-module franklab_msdrift `ml-config package_directory`/franklab_msdrift
ml-link-python-module franklab_mstaggedcuration `ml-config package_directory`/franklab_mstaggedcuration

### If no .mnt directory for day then run:
make_mda_ntrodeEpoch_links()

### If no .mountain directory for day then run nodejs script
node ./MS4setup_NTlinks.node.js '/data2/edeno/remy/preprocessing'

Set temp directory to copied data directory
Run Spike sorting
