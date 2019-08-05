import errno
import glob
import json
import logging
import os
import subprocess

import franklab_mountainsort.ms4_franklab_pyplines as pyp

logging.basicConfig(level='INFO', format='%(asctime)s %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S')


def move_mda_data(source_animal_path, target_animal_path, animal, dates):
    '''Move data from preprocessing to scratch.

    Parameters
    ----------
    source_animal_path : str
    target_animal_path : str
    animal : str
    dates : list of int

    '''
    for date in dates:
        date = str(date)
        source_mda_paths = glob.glob(
            os.path.join(source_animal_path, date, '*.mda'))
        source_mda_paths.sort()
        for source_path in source_mda_paths:
            target_dir = os.path.join(target_animal_path, date)
            os.makedirs(target_dir, exist_ok=True)
            logging.info(f'\nCopying {source_path} to {target_dir}\n')
            subprocess.run(
                f'rsync -avP {source_path} {target_dir}', shell=True)


def run_spike_sorting(animal, dates, ntrodes, input_path, output_path,
                      metrics_input='metrics_merged.json',
                      metrics_output='metrics_merged_tagged.json',
                      firing_rate_thresh=0.01, isolation_thresh=0.95,
                      noise_overlap_thresh=0.03, peak_snr_thresh=1.5,
                      extract_marks=True, extract_clips=True,
                      clip_size=100, freq_min=300, freq_max=6000,
                      adjacency_radius=-1, detect_threshold=3, detect_sign=-1,
                      sampling_rate=30000):
    '''Runs mountain sort on data.

    Parameters
    ----------
    animal : str
    dates : list of str
    ntrodes : list of int
    input_path : str
    output_path : str
    metrics_input : str, optional
    metrics_output : str, optional
    firing_rate_thresh : float, optional
    isolation_thresh : float, optional
    noise_overlap_thresh : float, optional
    peak_snr_thresh : float, optional
    extract_marks : bool, optional
    extract_clips : bool, optional
    clip_size : float, optional
         The size of extract clips around each spike in samples.
    freq_min : float, optional
        The highpass or low cutoff of the filter in Hz.
    freq_max : float, optional
        The lowpass or high cutoff of the filter in Hz.
    adjacency_radius : float, optional
        The radius in µm that defines a neighborhood of channels on which to
        sort (default -1 to ignore and not require a geom.csv file, useful for
        tetrodes).
    detect_threshold : float, optional
        Spike detection threshold in standard deviations.
    detect_sign : int, optional
         The direction of spike to detect (-1 for negative, 1 for positive,
         0 for both). -1 is recommended for most recordings.
    sampling_rate : int, optional
        Number of samples per second.

    '''
    # Make .mnt symlinks if they don't exist
    date_folders = os.path.join(input_path, '*')
    if len(glob.glob(os.path.join(date_folders, '*.mnt'))) == 0:
        make_mda_ntrodeEpoch_links(input_path)

    for date in dates:
        date = str(date)
        logging.info(f'Running {animal} date: {date} ntrodes: {ntrodes}')
        mountain_mda_path = os.path.join(
            input_path, date, f'{date}_{animal}.mountain')
        mountain_out_path = os.path.join(output_path, date, 'ms4')

        for electrode_number in ntrodes:
            mountain_mda_nt_path = os.path.join(
                mountain_mda_path, f'nt{electrode_number}')
            mountain_out_nt_path = os.path.join(
                mountain_out_path, f'nt{electrode_number}')
            os.makedirs(mountain_out_nt_path, exist_ok=True)
            mda_opts = {'anim': animal,
                        'date': date,
                        'ntrode': electrode_number,
                        'data_location': input_path}
            raw_mda = os.path.join(mountain_mda_nt_path, 'raw.mda')

            # create concatenated mda if it doesn't exist
            if not os.path.isfile(raw_mda):
                os.makedirs(mountain_mda_nt_path, exist_ok=True)
                # create params if it doesn't exist
                params_file = os.path.join(mountain_mda_nt_path, 'params.json')
                if not os.path.isfile(params_file):
                    params = {'samplerate': sampling_rate}
                    with open(params_file, 'w') as f:
                        json.dump(params, f, indent=4, sort_keys=True)
                logging.info(f'Creating concatenated epochs .mda: {raw_mda}')
                pyp.concat_epochs(dataset_dir=mountain_mda_nt_path,
                                  mda_opts=mda_opts)

            logging.info(f'####### NTRODE_INPUT: {raw_mda}')
            logging.info(f'####### NTRODE_OUTPUT: {mountain_out_nt_path}')

            pyp.filt_mask_whiten(
                dataset_dir=mountain_mda_nt_path,
                output_dir=mountain_out_nt_path,
                freq_min=freq_min,
                freq_max=freq_max)
            pyp.ms4_sort_on_segs(
                dataset_dir=mountain_mda_nt_path,
                output_dir=mountain_out_nt_path,
                adjacency_radius=adjacency_radius,
                detect_threshold=detect_threshold,
                detect_sign=detect_sign,
                mda_opts=mda_opts)
            pyp.merge_burst_parents(
                dataset_dir=mountain_mda_nt_path,
                output_dir=mountain_out_nt_path)
            pyp.add_curation_tags(
                dataset_dir=mountain_out_nt_path,
                output_dir=mountain_out_nt_path,
                metrics_input=metrics_input,
                metrics_output=metrics_output,
                firing_rate_thresh=firing_rate_thresh,
                isolation_thresh=isolation_thresh,
                noise_overlap_thresh=noise_overlap_thresh,
                peak_snr_thresh=peak_snr_thresh)

            if extract_marks:
                pyp.extract_marks(dataset_dir=mountain_out_nt_path,
                                  output_dir=mountain_out_nt_path)

            if extract_clips:
                pyp.extract_clips(dataset_dir=mountain_out_nt_path,
                                  output_dir=mountain_out_nt_path,
                                  clip_size=clip_size)


def make_mda_ntrodeEpoch_links(preprocessing_path):
    '''This makes the .mnt directory and its contents, containing symbolic
    links to the mda files for each epoch for each ntrode.

    We use this because our raw data is stored by epoch, with all tets in
    subfolders. Since we generally want to sort across all the epochs in a day,
    this tool reorganizes the data so that it is first grouped by tetrode, then
    epoch. This is all done with symlinks, so no duplicate files are generated.

    Parameters
    ----------
    preprocessing_path : str
        Folder path to animal preprocessing folder e.g.
        '.../{animal}/preprocessing'

    '''
    mda_filepaths = glob.glob(
        os.path.join(preprocessing_path, '*', '*', '*.nt*.mda'))
    for mda_filepath in mda_filepaths:
        try:
            mda_filename = os.path.basename(mda_filepath)
            date, animal, epoch, ending = mda_filename.split('_')
            _, electrode_name, _ = ending.split('.')
            dest_path = os.path.join(
                preprocessing_path, date, f'{date}_{animal}.mnt',
                f'{date}_{animal}.{electrode_name}.mnt')
            os.makedirs(dest_path, exist_ok=True)
            dest_link = os.path.join(dest_path, mda_filename)

            # to overwrite. remove ntlink if it already exists
            remove_NT_file(dest_link)

            # create directory of sym links to original mda
            os.symlink(mda_filepath, dest_link)

        except ValueError:
            pass


def remove_NT_file(ntrode_filename):
    try:
        os.remove(ntrode_filename)
    except OSError as e:
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occurred
