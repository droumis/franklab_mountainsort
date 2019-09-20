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
2. Create a DataFrame with .mda file information (`get_mda_files_dataframe`)
3. Run spike sorting on all or a subset of the DataFrame (`spike_sort_all`)
