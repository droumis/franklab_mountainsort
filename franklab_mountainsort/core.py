import glob
import json
import logging
import os
import subprocess

import dask
import pandas as pd

import franklab_mountainsort.ms4_franklab_pyplines as pyp


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


def spike_sort_all(mda_file_info, input_path, output_path,
                   metrics_input='metrics_merged.json',
                   metrics_output='metrics_merged_tagged.json',
                   firing_rate_thresh=0.01, isolation_thresh=0.95,
                   noise_overlap_thresh=0.03, peak_snr_thresh=1.5,
                   extract_marks=True, extract_clips=True,
                   clip_size=100, freq_min=300, freq_max=6000,
                   adjacency_radius=-1, detect_threshold=3, detect_sign=-1,
                   sampling_rate=30000):
    '''Runs mountain sort on all electrodes in `mda_file_info`

    Parameters
    ----------
    mda_file_info : pandas.DataFrame
        Dataframe generated using function `get_mda_files_dataframe`
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
    electrodes = mda_file_info.groupby(
        ['animal', 'date', 'electrode_number'])
    results = []
    for (animal, date, electrode_number), electrodes_df in electrodes:
        results.append(
            spike_sort_electrode(
                animal, date, electrode_number, input_path,
                output_path, metrics_input, metrics_output,
                firing_rate_thresh, isolation_thresh,
                noise_overlap_thresh, peak_snr_thresh,
                extract_marks, extract_clips, clip_size, freq_min,
                freq_max, adjacency_radius, detect_threshold,
                detect_sign, sampling_rate))
    dask.compute(*results)


@dask.delayed
def spike_sort_electrode(animal, date, electrode_number, input_path,
                         output_path, metrics_input='metrics_merged.json',
                         metrics_output='metrics_merged_tagged.json',
                         firing_rate_thresh=0.01, isolation_thresh=0.95,
                         noise_overlap_thresh=0.03, peak_snr_thresh=1.5,
                         extract_marks=True, extract_clips=True,
                         clip_size=100, freq_min=300, freq_max=6000,
                         adjacency_radius=-1, detect_threshold=3,
                         detect_sign=-1, sampling_rate=30000, geom=None):
    '''Runs mountain sort on all electrodes in `mda_file_info`

    Parameters
    ----------
    animal : str
    date : str
    electrode_number : int
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
    geom : ndarray or None, shape (n_contacts, 2), optional
        Geometry of the electrode. Important for probes.

    '''
    date = str(date)

    # Setup log file
    log_directory = os.path.join(output_path, 'logs')
    os.makedirs(log_directory, exist_ok=True)
    log_file = os.path.join(
        log_directory, f'{animal}_{date}_nt{electrode_number}.log')

    logger = logging.getLogger(f'{animal}_{date}_nt{electrode_number}')
    fhandler = logging.FileHandler(filename=log_file, mode='a')
    formatter = logging.Formatter(fmt='%(asctime)s %(message)s',
                                  datefmt='%d-%b-%y %H:%M:%S')
    fhandler.setFormatter(formatter)
    logger.addHandler(fhandler)
    logger.setLevel(logging.DEBUG)

    logger.info(
        f'Processing animal: {animal}, date: {date}, '
        f'electrode: {electrode_number}')
    logger.info(f'Parameters: {locals()}')

    mountain_mda_electrode_dir = os.path.join(
        input_path, date, f'{date}_{animal}.mountain',
        f'nt{electrode_number}')
    mda_opts = {'anim': animal,
                'date': date,
                'ntrode': electrode_number,
                'data_location': input_path}
    raw_mda = os.path.join(mountain_mda_electrode_dir, 'raw.mda')

    # Concatenate mda files from the same epoch
    if not os.path.isfile(raw_mda):
        logger.info('Concatenated .mda file does not exist.'
                    f'Creating {raw_mda}')
        os.makedirs(mountain_mda_electrode_dir, exist_ok=True)
        # create params if it doesn't exist
        params_file = os.path.join(
            mountain_mda_electrode_dir, 'params.json')
        if not os.path.isfile(params_file):
            params = {'samplerate': sampling_rate}
            with open(params_file, 'w') as f:
                json.dump(params, f, indent=4, sort_keys=True)
        pyp.concat_epochs(dataset_dir=mountain_mda_electrode_dir,
                          mda_opts=mda_opts)

    # Make sure .mountain output directory exists
    mountain_out_electrode_dir = os.path.join(
        output_path, date, 'ms4', f'nt{electrode_number}')
    os.makedirs(mountain_out_electrode_dir, exist_ok=True)

    logger.info('Preprocessing waveforms...')
    pyp.filt_mask_whiten(
        dataset_dir=mountain_mda_electrode_dir,
        output_dir=mountain_out_electrode_dir,
        freq_min=freq_min,
        freq_max=freq_max)

    logger.info('Sorting spikes...')
    pyp.ms4_sort_on_segs(
        dataset_dir=mountain_mda_electrode_dir,
        output_dir=mountain_out_electrode_dir,
        geom=geom,
        adjacency_radius=adjacency_radius,
        detect_threshold=detect_threshold,
        detect_sign=detect_sign,
        mda_opts=mda_opts)

    logger.info('Merging burst parents...')
    pyp.merge_burst_parents(
        dataset_dir=mountain_mda_electrode_dir,
        output_dir=mountain_out_electrode_dir)

    logger.info('Adding curation tags...')
    pyp.add_curation_tags(
        dataset_dir=mountain_out_electrode_dir,
        output_dir=mountain_out_electrode_dir,
        metrics_input=metrics_input,
        metrics_output=metrics_output,
        firing_rate_thresh=firing_rate_thresh,
        isolation_thresh=isolation_thresh,
        noise_overlap_thresh=noise_overlap_thresh,
        peak_snr_thresh=peak_snr_thresh)

    if extract_marks:
        logger.info('Extracting marks...')
        pyp.extract_marks(dataset_dir=mountain_out_electrode_dir,
                          output_dir=mountain_out_electrode_dir)

    if extract_clips:
        logger.info('Extracting clips...')
        pyp.extract_clips(dataset_dir=mountain_out_electrode_dir,
                          output_dir=mountain_out_electrode_dir,
                          clip_size=clip_size)
    logger.info('Done')


def get_mda_files_dataframe(data_path, recursive=False):
    '''

    Parameters
    ----------
    data_path : str
    recursive : bool, optional
        Recursive search using glob.

    Returns
    -------
    mda_files_dataframe : pandas.DataFrame

    Examples
    --------
    ```
    all_animal_info = get_mda_files_dataframe(
        '/data2/data1_backup/anna/*/preprocessing/')
    ```

    ```
    single_animal_info = get_mda_files_dataframe(
        '/data2/data1_backup/anna/remy/preprocessing/')
    ```
    '''

    mda_files = glob.glob(os.path.join(data_path, '*', '*.mda', '*.*.mda'),
                          recursive=recursive)
    file_info = [_get_file_information(mda_file) for mda_file in mda_files
                 if _get_file_information(mda_file) is not None]
    COLUMNS = ['animal', 'date', 'epoch', 'electrode_number', 'task',
               'relative_filepath']
    return (pd.DataFrame(file_info, columns=COLUMNS)
            .set_index(['animal', 'date', 'epoch', 'electrode_number'])
            .sort_index())


def _get_file_information(mda_file):
    try:
        date, animal, epoch, other = os.path.basename(mda_file).split('_')
        date, epoch = int(date), int(epoch)
        task, electrode_name, _ = other.split('.')
        electrode_number = int(electrode_name.strip('nt'))
        relative_filepath = os.path.join(
            animal, 'preprocessing', str(date), os.path.basename(mda_file))

        return animal, date, epoch, electrode_number, task, relative_filepath
    except ValueError:
        pass
