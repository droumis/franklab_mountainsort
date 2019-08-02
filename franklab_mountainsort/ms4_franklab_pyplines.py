#from mountainlab_pytools import mdaio
#from mountainlab_pytools import mlproc as mlp
import os
import json
import subprocess
import ms4_franklab_proc2py as p2p
import math
import re
from collections import defaultdict
# This script calls the helper functions defined in p2p that in turn, call MS processors
# This should be a collection of common processing steps that are standard across the lab, altho params can be changed flexibly
# Default params are defined in the arguments of each pypline, but can be overwritten by the user

# These pyplines should be called by a python batch script, which will manage running the steps on particular animal, days, and ntrodes
# AKGillespie based on code from JMagland
# Demetris Roumis



#before anything else, must concat all eps together becuase ms4 no longer handles the prv list of mdas

def concat_eps(*,dataset_dir, mda_list=[], opts={}, mda_opts={}):
    # runs 'ms3.concat_timeseries' using either
    # 1: mda_list provided
    # 2: mda_list empty and date, ntrode specified in opts
    # 3: string path to prv file containing entries for mda files

    strstart = []
    if type(mda_list) == list and mda_list != []:
        print('using provided list of mda files')
        for entry in mda_list:
            strstart.append('timeseries_list:'+entry)

    if mda_list == [] and set(['anim', 'date', 'ntrode', 'data_location']).issubset(mda_opts):
        print(f'scavenging list of mda file from mda directories of date:{mda_opts["date"]} ntrode:{mda_opts["ntrode"]}')
        mda_list = p2p.get_mda_list(mda_opts['anim'], mda_opts['date'], mda_opts['ntrode'], mda_opts['data_location'])

        for entry in mda_list:
            strstart.append('timeseries_list:'+entry)


    if type(mda_list) == str:
        print('using mda files listed in prv file')
        with open(mda_list) as f:
            mdalist=json.load(f)
        for entries in mdalist['files']:
            strstart.append('timeseries_list:'+entries['prv']['original_path'])

    joined = ' '.join(strstart)
    outpath = 'timeseries_out:'+dataset_dir+'/raw.mda'
#     print((['ml-run-process','ms3.concat_timeseries','--inputs',joined,'--outputs',outpath]))
    subprocess.call(['ml-run-process','ms3.concat_timeseries','--inputs',joined,'--outputs',outpath])

def filt_mask_whiten(*,dataset_dir,output_dir,freq_min=300,freq_max=6000,mask_artifacts=1,opts={}):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Dataset parameters
    ds_params=p2p.read_dataset_params(dataset_dir)

    # Bandpass filter
    p2p.bandpass_filter(
        timeseries=dataset_dir+'/raw.mda',
        timeseries_out=output_dir+'/filt.mda.prv',
        samplerate=ds_params['samplerate'],
        freq_min=freq_min,
        freq_max=freq_max,
        opts=opts
    )
    # Mask out artifacts
    if mask_artifacts:
        p2p.mask_out_artifacts(
            timeseries=output_dir+'/filt.mda.prv',
            timeseries_out=output_dir+'/filt.mda.prv',
            threshold = 5,
            interval_size=2000,
            opts=opts
            )
    # Whiten
    p2p.whiten(
        timeseries=output_dir+'/filt.mda.prv',
        timeseries_out=output_dir+'/pre.mda.prv',
        opts=opts
    )


# full = sort the entire file as one mda
def ms4_sort_full(*,dataset_dir, output_dir, geom=[], adjacency_radius=-1,detect_threshold=3,detect_sign=0,opts={}):

    # Fetch dataset parameters
    ds_params=p2p.read_dataset_params(dataset_dir)


    p2p.ms4alg(
        timeseries=output_dir+'/pre.mda.prv',
        geom=geom,
        firings_out=output_dir+'/firings_raw.mda',
        adjacency_radius=adjacency_radius,
        detect_sign=detect_sign,
        detect_threshold=detect_threshold,
        opts=opts
    )

    # Compute cluster metrics
    p2p.compute_cluster_metrics(
        timeseries=output_dir+'/pre.mda.prv',
        firings=output_dir+'/firings_raw.mda',
        metrics_out=output_dir+'/metrics_raw.json',
        samplerate=ds_params['samplerate'],
        opts=opts
    )


# segs = sort by timesegments, then join any matching  clusters
def ms4_sort_on_segs(*,dataset_dir, output_dir, geom=[], adjacency_radius=-1,detect_threshold=3,detect_sign=0,rm_segment_intermediates=1, opts={}, mda_opts={}):

    # Fetch dataset parameters
    ds_params=p2p.read_dataset_params(dataset_dir)

    if set(['anim', 'date', 'ntrode', 'data_location']).issubset(mda_opts):
        print(f'scavenging list of mda file from mda directories of date:{mda_opts["date"]} ntrode:{mda_opts["ntrode"]}')
        mda_list = p2p.get_mda_list(mda_opts['anim'], mda_opts['date'], mda_opts['ntrode'], mda_opts['data_location'])
        # calculate time_offsets and total_duration
        sample_offsets, total_samples = p2p.get_epoch_offsets(dataset_dir=dataset_dir, opts={'mda_list':mda_list})

    else:
        # calculate time_offsets and total_duration
        sample_offsets, total_samples = p2p.get_epoch_offsets(dataset_dir=dataset_dir)

    #break up preprocesed data into segments and sort each
    firings_list=[]
    timeseries_list=[]
    for segind in range(len(sample_offsets)):
        t1=math.floor(sample_offsets[segind])
        if segind==len(sample_offsets)-1:
            t2=total_samples-1
        else:
            t2=math.floor(sample_offsets[segind+1])-1

        segment_duration = t2-t1
        print('Segment '+str(segind+1)+': t1='+str(t1)+', t2='+str(t2)+', t1_min='+str(t1/ds_params['samplerate']/60)+', t2_min='+str(t2/ds_params['samplerate']/60));

        pre_outpath= dataset_dir+'/pre-'+str(segind+1)+'.mda'
        p2p.pyms_extract_segment(
            timeseries=output_dir+'/pre.mda.prv',
            timeseries_out=pre_outpath,
            t1=t1,
            t2=t2,
            opts=opts)

        firings_outpath=dataset_dir+'/firings-'+str(segind+1)+'.mda'
        p2p.ms4alg(
            timeseries=pre_outpath,
            firings_out=firings_outpath,
            geom=geom,
            detect_sign=detect_sign,
            adjacency_radius=adjacency_radius,
            detect_threshold=detect_threshold,
            opts=opts)

        firings_list.append(firings_outpath)
        timeseries_list.append(pre_outpath)

    firings_out_final=output_dir+'/firings_raw.mda'
    # sample_offsets have to be converted into a string to be properly passed into the processor
    str_sample_offsets=','.join(map(str,sample_offsets))
    print(str_sample_offsets)

    p2p.pyms_anneal_segs(
        timeseries_list=timeseries_list,
        firings_list=firings_list,
        firings_out=firings_out_final,
        dmatrix_out=[],
        k1_dmatrix_out=[],
        k2_dmatrix_out=[],
        dmatrix_templates_out=[],
        time_offsets=str_sample_offsets
    )

    # clear the temp pre and firings files if specified
    if rm_segment_intermediates:
        p2p.clear_seg_files(
            timeseries_list=timeseries_list,
            firings_list=firings_list
        )

    # Compute cluster metrics
    p2p.compute_cluster_metrics(
        timeseries=output_dir+'/pre.mda.prv',
        firings=output_dir+'/firings_raw.mda',
        metrics_out=output_dir+'/metrics_raw.json',
        samplerate=ds_params['samplerate'],
        opts=opts
    )

def merge_burst_parents(*, dataset_dir, output_dir,opts={}):

    p2p.merge_burst_parents(
        firings=output_dir+'/firings_raw.mda',
        metrics=output_dir+'/metrics_raw.json',
        firings_out=output_dir+'/firings_burst_merged.mda',
        opts={}
    )

    ds_params=p2p.read_dataset_params(dataset_dir)
    # Compute cluster metrics
    p2p.compute_cluster_metrics(
        timeseries=output_dir+'/pre.mda.prv',
        firings=output_dir+'/firings_burst_merged.mda',
        metrics_out=output_dir+'/metrics_merged.json',
        samplerate=ds_params['samplerate'],
        opts={}
    )

def add_curation_tags(*, dataset_dir, output_dir, firing_rate_thresh=.01, isolation_thresh=.95, noise_overlap_thresh=.03, peak_snr_thresh=1.5, metrics_input = '', metrics_output = '', opts={}):
    # note that this is split out and not included after metrics calculation
    # because of a bug in ms3.combine_cluster_metrics - doesn't work if anything follows it
    if not metrics_input:
        metrics_input = '/metrics_raw.json'
    if not metrics_output:
        metrics_output = '/metrics_tagged.json'

    p2p.tagged_curation(
        cluster_metrics=dataset_dir+metrics_input,
        metrics_tagged=output_dir+metrics_output,
        firing_rate_thresh=firing_rate_thresh,
        isolation_thresh=isolation_thresh,
        noise_overlap_thresh=noise_overlap_thresh,
        peak_snr_thresh=peak_snr_thresh,
        mv2file='',
        opts=opts
    )


def recalc_metrics(*, dataset_dir, output_dir, firings_in= '', metrics_to_update = '', firing_rate_thresh=.01, isolation_thresh=.95, noise_overlap_thresh=.03, peak_snr_thresh=1.5, mv2_file='', opts={}):
    #post-merge, should recalculate metrics and update tags (both tags based on thresholds and any manually added ones, stored in the mv2)

    # untested!
    if not firings_in:
        firings_in = '/firings_processed.json'
    if not metrics_to_update:
        metrics_to_update = '/metrics_tagged.json'

#     ds_params=p2p.read_dataset_params(dataset_dir)

    p2p.compute_cluster_metrics(
        timeseries=output_dir+'/pre.mda.prv',
        firings=output_dir+firings_in,
        metrics_to_update = output_dir+metrics_to_update,
        samplerate=ds_params['samplerate'],
        opts={}
    )

    p2p.tagged_curation(
        cluster_metrics=dataset_dir+metrics_to_update,
        metrics_tagged=output_dir+metrics_to_update,
        firing_rate_thresh=firing_rate_thresh,
        isolation_thresh=isolation_thresh,
        noise_overlap_thresh=noise_overlap_thresh,
        peak_snr_thresh=peak_snr_thresh,
        mv2file=mv2_file,
        opts=opts
    )

def extract_clips(*,dataset_dir, output_dir, clip_size=100, opts={}):
    opts['clip_size'] = clip_size

    p2p.pyms_extract_clips(
        timeseries=dataset_dir+'/pre.mda.prv',
        firings=dataset_dir+'/firings_raw.mda',
        clips_out=output_dir+'/clips.mda',
        opts=opts)

def extract_marks(*,dataset_dir, output_dir, opts={}):
    opts['clip_size'] = 1

    p2p.pyms_extract_clips(
        timeseries=dataset_dir+'/pre.mda.prv',
        firings=dataset_dir+'/firings_raw.mda',
        clips_out=output_dir+'/marks.mda',
        opts=opts)

# def extract_marks(*, dataset_dir, output_dir, opts={}):
#     p2p.pyms_extract_marks(
#         timeseries=dataset_dir+'/pre.mda.prv',
#         firings=dataset_dir+'/firings_raw.mda',
#         marks_out=output_dir+'/marks.mda',
#         markstimes_out=output_dir+'/markstimes.mda',
#         opts=opts)
