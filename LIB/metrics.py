#!/usr/bin/env python
#import sys
#sys.path.append('/Users/thompsong/src/kitchensinkGT/LIB')
import os
import numpy as np
from scipy.stats import describe
from obspy import Stream
from obspy.signal.quality_control import MSEEDMetadata 
#import libseisGT
from libseisGT import add_to_trace_history, clean_trace
import pandas as pd

"""
Functions for computing data quality metrics and statistical metrics (such as amplitude, energy and frequency) 
on Stream/Trace objects.

In terms of order of application:

1. Read the raw data.
2. Fix Trace IDs.
3. Compute QC metrics (and potentially remove bad Traces).
4. Correct the data (and save corrected data as MSEED/StationXML).
5. Compute statistical metrics.

"""


def process_trace(tr, inv=None):
# function tr, quality_factor, snr = compute_metrics(tr)
# This function wraps all others
    # Here we compute simple metrics on each trace (and used to write them to NET.STA.CHAN files). 
    # These metrics are:
    #     1. duration of signal
    #     2. signal amplitude
    #     3. noise amplitude
    #     4. signal-to-noise ratio
    
    if not 'history' in tr.stats:
        tr.stats['history'] = list()     
        
    """ RAW DATA QC METRICS """
    qcTrace(tr)
    tr.stats["quality_factor"] = trace_quality_factor(tr) #0 = blank trace, 1 = has some 0s and -1s, 3 = all looks good
    tr.stats.quality_factor -= tr.stats.metrics['num_gaps']
    tr.stats.quality_factor -= tr.stats.metrics['num_overlaps']
    tr.stats.quality_factor -= (100.0 - tr.stats.metrics['percent_availability'] )

    if tr.stats.quality_factor <= 0:
        return

    """ CLEAN (DETREND, BANDPASS, CORRECT) TRACE """
    clean_trace(tr, taperFraction=0.05, filterType="bandpass", freq=[0.1, 40.0], corners=2, zerophase=True, inv=inv)

    """ CLEAN DATA WAVEFORM METRICS """ 
    # estimate signal-to-noise ratio (after detrending)
    tr.stats["snr"] = signaltonoise(tr)
    add_to_trace_history(tr, 'Signal to noise measured.')
    tr.stats["quality_factor"] += np.log10(tr.stats["snr"][0])
    
    tr.stats["duration"] = tr.stats.npts /  tr.stats.sampling_rate # before or after detrending

    # SciPy stats
    #tr.stats["scipy"] = describe(tr.data, nan_policy = 'omit') # only do this after removing trend/high pass filtering
    tr.stats["scipy"] = describe(tr.data, nan_policy = 'omit')._asdict()
    add_to_trace_history(tr, 'scipy.stats metrics added as tr.stats.scipy.')
    
    # Update other stats
    qcTrace(tr)
    
    
def qcTrace(tr):
    """ qcTrace(tr) DATA QUALITY CHECKS """
    
    """ Useful MSEED metrics
    {'start_gap': None, 'end_gap': None, 'num_gaps': 0, 
     'sum_gaps': 0, 'max_gap': None, 'num_overlaps': 0, 
     'sum_overlaps': 0, 'max_overlap': None, 'quality': 'D', 
     'sample_min': -22404, 'sample_max': 9261, 
     'sample_mean': -3854.7406382978725, 'sample_median': -3836.0, 
     'sample_lower_quartile': -4526.0, 'sample_upper_quartile': -3105.0, 
     'sample_rms': 4426.1431329789848, 
     'sample_stdev': 2175.2511682727431, 
     'percent_availability': 100.0}           
    """
    tr.write('temp.mseed')
    mseedqc = MSEEDMetadata(['temp.mseed']) 
    tr.stats['metrics'] = mseedqc.meta
    os.remove('temp.mseed')
    add_to_trace_history(tr, 'MSEED metrics computed (similar to ISPAQ/MUSTANG).')

def _detectClipping(tr, countThresh = 10):
    upper_clipped = False
    lower_clipped = False
    y = tr.data
    mu = np.nanmax(y)
    md = np.nanmin(y)
    countu = (tr.data == mu).sum()
    countd = (tr.data == md).sum()
    if countu >= countThresh:
        add_to_trace_history(tr, 'Trace %s appears to be clipped at upper limit %e (count=%d)' % (tr.id, mu, countu) )    
        upper_clipped = True
    if countd >= countThresh:
        add_to_trace_history(tr, 'Trace %s appears to be clipped at lower limit %e (count=%d)' % (tr.id, mu, countu) )       
        lower_clipped = True
    return upper_clipped, lower_clipped
    
def trace_quality_factor(tr):
    # trace_quality_factor(tr)
    # a good trace has quality factor 3, one with 0s and -1s has 1, bad trace has 0
    quality_factor = 1
    
    # ignore traces with few samples
    if tr.stats.npts < 100:
        add_to_trace_history(tr, 'Not enough samples')
        return 0
    
    # ignore traces with weirdly low sampling rates
    if tr.stats.sampling_rate < 19.99:
        add_to_trace_history(tr, 'Sampling rate too low')
        return 0

    # ignore blank trace
    anyData = np.count_nonzero(tr.data)
    if anyData==0:
        add_to_trace_history(tr, 'Trace is blank')
        return 0
    
    # check for bit level noise
    u = np.unique(tr.data)
    num_unique_values = u.size
    if num_unique_values > 10:
        quality_factor += np.log10(num_unique_values)
    else:
        add_to_trace_history(tr, 'bit level noise suspected')
        return 0

    # check for sequences of 0 or 1
    trace_good_flag = _check0andMinus1(tr.data)
    if trace_good_flag:
        quality_factor += 1
    else:
        add_to_trace_history(tr, 'sequences of 0 or -1 found')
        
    # check if trace clipped
    upperClipped, lowerClipped = _detectClipping(tr) # can add another function to interpolate clipped values
    if not upperClipped:
        quality_factor += 1
    if not lowerClipped:
        quality_factor += 1
         
    # check for outliers
    outlier_count, outlier_indices = _mad_based_outlier(tr, thresh=50.0)
    print('Outliers: %d' % outlier_count)
    if outlier_count == 0:
        quality_factor += 1
    else:
        add_to_trace_history(tr, '%d outliers found' % outlier_count)
        tr.stats['outlier_indices'] = outlier_indices    
       
    return quality_factor    

def _mad_based_outlier(tr, thresh=3.5):
    tr2 = tr.copy()
    tr2.detrend()
    points = tr2.data
    if len(points.shape) == 1:
        points = points[:,None]
    #points = np.absolute(points)
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation
    
    outlier_indices = np.array(np.where(modified_z_score > thresh))
    outlier_count = outlier_indices.size
    if outlier_count > 0:
        print('size diff = %d, median = %e, med_abs_deviation = %e ' % (diff.size, median, med_abs_deviation))
        mzs = sorted(modified_z_score)
        print(mzs[-10:])
    
    return outlier_count, outlier_indices

def _check0andMinus1(liste):
# function bool_good_trace = check0andMinus1(tr.data)
    liste=list(liste)
    listStr=''.join(str(i) for i in liste)
    if  "000000000000" in listStr or "-1-1-1-1-1-1-1-1" in listStr :
        return False
    else:
        return True   
    
def signaltonoise(tr):
# function snr, highval, lowval = signaltonoise(tr)
    # Here we just make an estimate of the signal-to-noise ratio
    #
    # Normally the trace should be pre-processed before passing to this routine, e.g.
    # * remove ridiculously large values
    # * remove any sequences of 0 from start or end
    # * detrend
    # * bandpass filter
    #
    # Processing:
    #    1. ensure we still have at least 10 seconds
    #    2. take absolute values
    #    3. compute the maximum of each 1-s of data, call this time series M
    #    4. compute 95th and 5th percentile of M, call these M95 and M5
    #    5. estimate signal-to-noise ratio as M95/M5
    
    highval = -1
    lowval = -1
    snr = -1 
    a = tr.data

    fsamp = int(tr.stats.sampling_rate)
    npts = tr.stats.npts
    numseconds = int(npts/fsamp)
    if numseconds > 10:
        a = a[0:int(fsamp * numseconds - 1)]             # remove any fractional second from end of trace
        abs = np.absolute(a)                             # take absolute value of a        
        abs_resize = np.resize(abs, (fsamp, numseconds)) # resize so that each second is one row
        M = np.max(abs_resize,axis=0)                    # find max of each second / row
        highval = np.nanpercentile(M,95)                    # set highval to 95th percentile, to represent signal amplitude
        if highval < 1:
            highval = -1
            return (snr, highval, lowval)
        lowval = np.nanpercentile(M,5)                      # set lowval to 5th percentile, to represent noise amplitude
        snr = highval / lowval
        print(abs_resize.shape)
        print(M.shape)
    return (snr, highval, lowval,)





def ssam(tr, f, S):
    if not f.size==S.size:
        return
    # use actual amplitude, not dB. 
    ssamValues = []
    tr.stats['spectral_amplitude'] = S
    tr.stats['spectral_frequencies'] = f
    for fmin in np.arange(0.0, 16.0, 1.0):
        f_indexes = np.intersect1d(np.where(f>=fmin), np.where(f<fmin+1.0))
        S_selected = S[f_indexes]
        ssamValues.append(np.nanmean(S_selected) )
    tr.stats['ssam'] = ssamValues 
    
def ampengfft(st):
    df = eventStatistics(st)
    #print(df)
    
    for tr in st:
        row = df[df['id'] == tr.id]
        tr.stats['peaktime'] = row.iloc[0]['time']
        tr.stats['peakamp'] = row.iloc[0]['peakamp']
        tr.stats['energy'] = row.iloc[0]['energy']
        
        # add SSAM
        if not 'ssam' in tr.stats:
            secsPerFFT = np.ceil((tr.stats.delta * tr.stats.npts)/100)
            [f, t, Zxx, A] = computeSingleSpectrogram(tr, secsPerFFT)
            try:
                S = np.mean(np.abs(Zxx), axis=0)
                ssam(tr, f, S)
            except:
                print(f.shape, Zxx.shape, A.shape, S.shape)
                return
            
def peak_amplitudes(st):   
    
    seismic1d_list = []
    seismic3d_list = []
    infrasound_list = []
    
    #ls.clean_stream(st, taper_fraction=0.05, freqmin=0.05, causal=True)
               
    # velocity, displacement, acceleration seismograms
    stV = st.select(channel='[ESBH]H?') 
    stD = stV.copy().integrate()
    for tr in stD:
        add_to_trace_history(tr, 'integrated')    
    stA = stV.copy().differentiate()
    for tr in stA:
        add_to_trace_history(tr, 'differentiated') 
     
    # Seismic vector data  
    stZ = stV.select(channel="[ESBH]HZ")
    for tr in stZ:
        thisID = tr.id[:-1]
        st3cV = stV.select(id = '%s[ENZ12RT]' % thisID)
        if len(st3cV)==3:
            st3cD = stD.select(id = '%s[ENZ12RT]' % thisID)
            st3cA = stA.select(id = '%s[ENZ12RT]' % thisID)
            md = ls.max_3c(st3cD)
            mv = ls.max_3c(st3cV)
            ma = ls.max_3c(st3cA)
            d = {'traceID':thisID, 'PGD':md[0], 'PGV':mv[0], 'PGA':ma[0], 'calib':tr.stats.calib, 'units':tr.stats.units}
            seismic3d_list.append(d)              
    seismic3d = pd.DataFrame(seismic3d_list)
    
    # Seismic 1-c data
    peakseismo1cfile = os.path.join(eventdir, 'summary_seismic_1c.csv')
    for c in range(len(stV)):
        md = max(abs(stD[c].data))
        mv = max(abs(stV[c].data))
        ma = max(abs(stA[c].data))  
        d = {'traceID':stV[c].id, 'PGD':md[0], 'PGV':mv[0], 'PGA':ma[0], 'calib':stV[c].stats.calib, 'units':stV[c].stats.units}
        seismic1d_list.append(d)    
    seismic1d = pd.DataFrame(seismic1d_list)        
            
    # Infrasound data
    peakinfrafile = os.path.join(eventdir, 'summary_infrasound.csv')
    stP = st.select(channel="[ESBH]D?")
    stPF = stP.copy().filter('bandpass', freqmin=1.0, freqmax=20.0, corners=2, zerophase=True)    
    for c in range(len(stP)):
        mb = max(abs(stP[c].data))
        mn = max(abs(stPF[c].data)) 
        d = {'traceID':stP[c].id, 'PP':mb[0], 'PPF':mn[0], 'calib':stP[c].stats.calib, 'units':stP[c].stats.units}
        infrasound_list.append(d)  
    infrasound = pd.DataFrame(infrasound_list)    
    
    return (seismic3d, seismic1d, infrasound)



def eventStatistics(st):
    """ 
    Take a Stream object, measure some basic metrics 
    
    Example: 
        eventStatistics(Stream)

    Returns: 
        list of dictionaries (one per trace)  with keywords id, sample, time, peakamp, energy


    Created for MiamiLakes project.
    
    Note: replace this with method from STA/LTA detection

    """
    # create an empty list of dictionaries
    list_of_dicts = []

    # for each trace, add a new dictionary to the list
    for tr in st:
        thisrow = dict()  # each dictionary will become a row of a dataframe

        # if we use absolute function, we don't need to check if maxy > -miny
        tr.detrend()
        y = np.absolute(tr.data)

        # we did this before
        maxy = y.max()
        maxy_i = y.argmax()
        maxy_t = tr.stats.starttime + maxy_i / tr.stats.sampling_rate

        # add new elements to dictionary
        thisrow['id'] = tr.id
        thisrow['sample'] = maxy_i
        thisrow['time'] = maxy_t
        thisrow['peakamp'] = maxy

        # add a new measurement: energy
        thisrow['energy'] = np.sum(np.square(y)) / tr.stats.sampling_rate

        # add row (dict) to list
        list_of_dicts.append(thisrow)

    # Convert list of dicts to dataframe
    df = pd.DataFrame.from_dict(list_of_dicts)

    return df

def choose_best_traces(st, MAX_TRACES=8, include_infrasound=False):
    # st2 = limit_traces(st, MAX_TRACES=8)

    if len(st)>MAX_TRACES:
        priority = [float(tr.stats.quality_factor) for tr in st]
        component = [tr.stats.channel[2] for tr in st]        
        for i, tr in enumerate(st):
            if component=="Z":
                priority[i] *= 2
            if tr.stats.channel[1]=='D':
                if include_infrasound:
                    priority[i] *= 2 
                else:
                    priority[i] *= 0                    
        
        j = np.argsort(priority)
        chosen = j[-MAX_TRACES:]
        
        #print(chosen)
        return chosen

def select_by_index_list(st, chosen):
    st2 = Stream()
    for i, tr in enumerate(st):
        if i in chosen:
            st2.append(tr)
    return st2