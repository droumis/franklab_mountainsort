
'''
debug mountainsort from script instead of ipynb
'''
print('hello world')

from franklab_mountainsort import get_mda_files_dataframe
from franklab_mountainsort import spike_sort_all
import logging
import os

## set animal
data_path = "/media/droumis/DR_swapdata9/"
animal = "JZ4"

## get mda file df
preprocessing_folder = os.path.join(data_path, animal, "preprocessing/")
print(f"data_path = {preprocessing_folder}")
mda_file_info = get_mda_files_dataframe(preprocessing_folder, recursive=True)

## set logging
FORMAT = "%(asctime)s %(message)s"
logging.basicConfig(level="INFO", format=FORMAT, datefmt="%d-%b-%y %H:%M:%S")

## set tmp dir
temp_path = os.path.join(data_path, animal, "mountainsort_tmp")
os.environ["ML_TEMPORARY_DIRECTORY"] = temp_path

##
spike_sort_all(
    mda_file_info.query("date == 20170424 and electrode_number == 28"),
    mountainlab_output_folder=None,
    firing_rate_thresh=0.01,
    isolation_thresh=0.95,
    noise_overlap_thresh=0.03,
    peak_snr_thresh=1.5,
    extract_marks=True,
    extract_clips=True,
    clip_time=1.5,
    freq_min=600,
    freq_max=6000,
    adjacency_radius=-1,
    detect_threshold=3,
    detect_interval=10,
    detect_sign=-1,
    sampling_rate=30000,
    drift_track=True,)