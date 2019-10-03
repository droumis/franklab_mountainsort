# franklab_mountainsort

### Installation
1. Check if miniconda or anaconda is installed
https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

2. Clone git repository
```bash
git clone https://github.com/edeno/franklab_mountainsort.git
```
3. Create conda environment
```bash
conda env create -f environment.yml
```
4. activate conda environment
```bash
conda activate franklab_mountainsort
```
5. Install package code
```bash
python setup.py install
```
6. Check if everything installed correctly. Open jupyter notebook, jupyter lab, or python console. Try to import franklab_mountainsort.
```python
import franklab_mountainsort
```

### Usage
1. Define input data directory, output data, temp directory
2. Create a pandas dataframe with `.mda` file information (`get_mda_files_dataframe`)
3. Run spike sorting on all or a subset of the DataFrame (`spike_sort_all`)
4. Open qt-mountainview (make sure you've activated the conda environment. If not type `conda activate franklab_mountainsort` into the terminal.)
```
qt-mountainview --pre=pre.mda.prv --firings=firings_raw.mda --samplerate=30000 --cluster_metrics=metrics_raw.json
```

+ raw.mda is in the `<animal>/preprocessing/<date>/<date>_<animal>.mountain`  folder
+ need to categorize all the files and where they are

### Inputs to MountainSort
+ `raw.mda.prv`: The concatenated time series for one day of recording.
+ `params.json`: contains information about the parameters to use for the sort. These can be overwritten by specifying them in the call to run the sort (in the batch script)
+ `geom.csv` (optional): contains information about the location of contacts for that ntrode; used in concert with adjacency_radius to determine the neighborhoods to sort on. In the case of tetrodes, this is not necessary because all the contacts of a tetrode should be sorted together as a single neighborhood. This can be specified by not providing a geom.csv, setting adjacency_radius to -1, or both.

### Note: `.prv` vs.`.mda` files
A `.prv` file is a text file with pointers to the binary file.
A `.mda` file is a binary file containing the actual data.

They can be used interchangeably when using `qt-mountainview`.


### Outputs from MountainSort
For each electrode:
+ `clips.mda`: The waveforms around the spike times.
+ `filt.mda.prv`: The time series after it has been bandpass filtered.
+ `firings_burst_merged.mda`: The firings file after burst merge processing.
+ `firings_raw.mda`: This contains the actual spike timing info that you care most about [electrode; time;label x #events….] for ALL detected events, regardless of cluster quality
+ `metrics_merged.json`: this contains the metrics for the curated clusters, such as isolation scores, noise overlap, SNR, and more.
+ `metrics_merged_tagged.json`: Same as metrics merged but with tags.
+ `metrics_raw.json`: metrics for all the original clusters
+ `pre.mda.prv`: The time series after it has been bandpass filtered and whitened.
