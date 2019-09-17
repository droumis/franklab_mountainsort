import os

from franklab_mountainsort import (get_mda_files_dataframe, move_mda_data,
                                   spike_sort_all)

import logging

FORMAT = '%(asctime)s %(message)s'

logging.basicConfig(level='INFO', format=FORMAT, datefmt='%d-%b-%y %H:%M:%S')

animal = 'kf19'
dates = ['20170913']
target_path = '/data2/edeno/'

source_animal_path = f'/data2/jason/{animal}/preprocessing'

input_path = os.path.join(target_path, animal, 'preprocessing')
output_path = os.path.join(target_path, animal, 'mountainlab_output')
temp_path = os.path.join(target_path, 'temp2')
log_directory = os.path.join(target_path, 'temp2')

os.environ['ML_TEMPORARY_DIRECTORY'] = temp_path

move_mda_data(
    source_animal_path, input_path, animal, dates)

mda_file_info = get_mda_files_dataframe(input_path)

spike_sort_all(
    mda_file_info, input_path, output_path,
    metrics_input='metrics_merged.json',
    metrics_output='metrics_merged_tagged.json',
    firing_rate_thresh=0.01, isolation_thresh=0.95,
    noise_overlap_thresh=0.03, peak_snr_thresh=3,
    extract_marks=False, freq_min=600)
