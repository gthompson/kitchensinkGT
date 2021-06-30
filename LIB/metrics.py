#!/usr/bin/env python
#import sys
#sys.path.append('/Users/thompsong/src/kitchensinkGT/LIB')
import libseisGT
import numpy as np

def check0andMinus1(liste):
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

def compute_metrics(tr):
# function tr, quality_factor, snr = compute_metrics(tr)
# This function wraps all others
    # Here we compute simple metrics on each trace and write them to NET.STA.CHAN files. 
    # These metrics are:
    #     1. duration of signal
    #     2. signal amplitude
    #     3. noise amplitude
    #     4. signal-to-noise ratio

    duration = tr.stats.npts /  tr.stats.sampling_rate
    quality_factor = trace_quality_factor(tr)
    snr = (-1, -1, -1)
    if quality_factor > 0:

        clean_trace(tr)
        correct_nslc(tr)

        # estimate signal-to-noise ratio
        snr = signaltonoise(tr)
        if snr[0] <= 1:
            quality_factor = quality_factor * 0.5

    return (tr, quality_factor, snr)

def trace_quality_factor(tr):
# function quality_factor = trace_quality_factor(tr)
    quality_factor = 1

    # ignore blank trace
    anyData = np.count_nonzero(tr.data)
    if anyData==0:
        quality_factor = 0 
        return quality_factor

    # check for sequences of 0 or 1
    trace_good_flag = check0andMinus1(tr.data)

    if not trace_good_flag:
        quality_factor = 0
        return quality_factor

    return quality_factor

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
            d = dict('traceID':thisID, 'PGD':md[0], 'PGV':mv[0], 'PGA':ma[0], 'calib':tr.stats.calib, 'units':tr.stats.units)
            seismic3d_list.append(d)              
    seismic3d = pd.DataFrame(seismic3d_list)
    
    # Seismic 1-c data
    peakseismo1cfile = os.path.join(eventdir, 'summary_seismic_1c.csv')
    for c in range(len(stV)):
        md = max(abs(stD[c].data))
        mv = max(abs(stV[c].data))
        ma = max(abs(stA[c].data))  
        d = dict('traceID':stV[c].id, 'PGD':md[0], 'PGV':mv[0], 'PGA':ma[0], 'calib':stV[c].stats.calib, 'units':stV[c].stats.units)
        seismic1d_list.append(d)    
    seismic1d = pd.DataFrame(seismic1d_list)        
            
    # Infrasound data
    peakinfrafile = os.path.join(eventdir, 'summary_infrasound.csv')
    stP = st.select(channel="[ESBH]D?")
    stPF = stP.copy().filter('bandpass', freqmin=1.0, freqmax=20.0, corners=2, zerophase=True)    
    for c in range(len(stP)):
        mb = max(abs(stP[c].data))
        mn = max(abs(stPF[c].data)) 
        d = dict('traceID':stP[c].id, 'PP':mb[0], 'PPF':mn[0], 'calib':stP[c].stats.calib, 'units':stP[c].stats.units)
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