import os

from franklab_mountainsort import (get_mda_files_dataframe, move_mda_data,
                                   spike_sort_all)

animal = 'remy'
dates = ['20170920', '20170921', '20170922']
target_path = '/data2/edeno/'

source_animal_path = f'/data2/data1_backup/anna/{animal}/preprocessing'

input_path = os.path.join(target_path, animal, 'preprocessing')
output_path = os.path.join(target_path, animal, 'mountainlab_output')
temp_path = os.path.join(target_path, 'temp')
log_directory = os.path.join(target_path, 'temp')

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
