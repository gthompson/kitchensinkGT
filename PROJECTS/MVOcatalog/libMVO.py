#!/usr/bin/env python

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

 
def clip_trace(tr):     
# function tr = clip_trace(tr)
    # remove absurdly large values
    AMP_LIMIT = 10000000
    a = tr.data
    np.clip(a, -AMP_LIMIT, AMP_LIMIT, out=a)
    np.where(a == AMP_LIMIT, 0, a)
    tr.data = a

def change_last_sample(tr):
# function tr = change_last_sample(tr)
    # For some SEISAN files from Montserrat - possibly from SEISLOG conversion,
    # the last value in the time series was always some absurdly large value
    # This sections was to change last value to be the mean of the rest
    #a = tr.data
    #a = np.delete(a,[np.size(a) - 1])
    #m = np.mean(a)
    #tr.data[-1] = m
    return tr

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

def swap32(i):
# function y = swap(tr.data)
    return struct.unpack("<i", struct.pack(">i", i))[0]

def fix_trace_id(st):
    for tr in st:
        nslc = correct_nslc(tr.id, tr.stats.sampling_rate)
        tr.id = nslc

def correct_nslc(traceID, Fs):
    # Montserrat trace IDs are often bad. return correct trace ID
    oldnet, oldsta, oldloc, oldcha = traceID.split('.')

    net = 'MV'
    chan = oldcha
    sta = oldsta
    loc = oldloc
    #sta = oldsta.strip()
    #chan = oldcha.strip()
    
    if loc == 'J':
        loc = ''
        
    if loc == '--':
        loc = ''
                
    # channel code is bandcode + instrumentcode + orientationcode
    bandcode = libseisGT.get_bandcode(Fs)
    
    instrumentcode = 'H' # seismic velocity sensor is default
        
    if chan[0]=='P' or chan[0:2]=='AP':
        instrumentcode ='D' # infrasound/acoustic
        orientationcode = 'F'
        if chan[1].isnumeric(): # e.g. channel like P5
            loc = chan[1]
        
    if instrumentcode == 'H': # seismic
        if len(chan)==2:
            orientationcode = loc
            loc = ''
        else:
            orientationcode = chan[2]

    chan = bandcode + instrumentcode + orientationcode

    return net + "." + sta + "." + loc + "." + chan

def detectClipping(tr):
    countThresh = 10
    y = tr.data
    mu = np.nanmax(y)
    md = np.nanmin(y)
    countu = (tr.data == mu).sum()
    countd = (tr.data == md).sum()
    if countu + countd >= countThresh:
        print('Trace %s appears to be clipped from %d (count=%d) to %d (count=%d)' % (tr.id, mu, countu, md, countd) )
        return True
    else:
        return False
        

def clean_trace(tr, taperFraction=0.05, filterType="bandpass", freq=[0.1, 20.0], corners=2, zerophase=True, inv=None):

    traceIsClipped = detectClipping(tr) # can add another function to interpolate clipped values
    
    # remove absurd values
    clip_trace(tr)
    
    # save the start and end times for later 
    startTime = tr.stats.starttime
    endTime = tr.stats.endtime
        
    # pad the Trace
    y = tr.data
    npts = tr.stats.npts
    npts_pad = int(taperFraction * npts)
    npts_pad_seconds = npts_pad * tr.stats.delta
    if npts_pad_seconds < 10.0: # impose a minimum pad length of 10-seconds
        npts_pad = int(10.0 / tr.stats.delta)
    y_prepend = np.flip(y[0:npts_pad])
    y_postpend = np.flip(y[-npts_pad:])
    y = np.concatenate( [y_prepend, y, y_postpend ] )
    padStartTime = startTime - npts_pad * tr.stats.delta
    tr.data = y
    tr.stats.starttime = padStartTime
    max_fraction = npts_pad / tr.stats.npts
    
    # clean
    tr.detrend('linear')
    tr.taper(max_percentage = max_fraction)
    
    if filterType == 'bandpass':
        tr.filter(filterType, freqmin=freq[0], freqmax=freq[1], corners=corners, zerophase=zerophase)
    else:    
        tr.filter(filterType, freq=freq, corners=corners, zerophase=zerophase) 
        
    # deconvolve
    if inv:
        if tr.stats.channel[1]=='H': # a seismic velocity channel
            tr.remove_response(inventory=inv, output="VEL")  
    
    # remove the pad
    tr.trim(starttime=startTime, endtime=endTime)

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
            
            
def clean_stream(st, inventory=None, remove_response_flag=False,
                   water_level=60, filter_flag=False, pre_filt=None,
                   starttime=None, endtime=None,
                   resample_flag=False, sampling_rate=100.0,
                   taper_type="hann", taper_percentage=0.05,
                   filterType="bandpass", freq=[0.1, 20.0], corners=2, zerophase=True ):
    """

    Based on code by Wenjie Lei (lei@princeton.edu), 2016. GNU Lesser General Public License, version 3 (LGPLv3)

    :param st: input stream
    :type st: obspy.Stream
    :param remove_response_flag: flag for remove instrument response. If True,
        then inv should be specified, and filter_flag would not be taken caren
        of. If you want just filter the seismogram, please leave this to False
        and set filter_flag to True.
    :type remove_response_flag: bool
    :param inventory: station inventory information
    :type inventory: obspy.Inventory
    :param water_level: water level used in remove instrument response. The
        default value in obspy is 60.
    :type water_level: float
    :param filter_flag:flag for filter the seismogram
    :type filter_flag: bool
    :param pre_filt: list of tuple of 4 corner frequency for filter,
        in ascending order(unit: Hz)
    :type pre_filt: list, tuple or numpy.array
    :param resample_flag: flag for data resampling
    :type resample_flag: bool
    :param sampling_rate: resampling rate(unit: Hz)
    :type sampling_rate: float
    :param taper_type: taper type, options from obspy taper
    :type taper_type: str
    :param taper_percentage: percentage of taper
    :type taper_percentage: float
    :return: processed stream
    """
    
    # check input data type
    if isinstance(st, Trace):
        st = Stream(traces=[st, ])
        _is_trace = True
    elif isinstance(st, Stream):
        _is_trace = False
    else:
        raise TypeError("Input seismogram should be either obspy.Stream "
                        "or obspy.Trace")

    if filter_flag or remove_response_flag:
        # detrend ,demean, taper
        st.detrend("linear")
        st.detrend("demean")
        st.taper(max_percentage=taper_percentage, type=taper_type)

    # remove response or filter
    if filter_flag:
        if pre_filt is None or len(pre_filt) != 4:
            raise ValueError("Filter band should be list or tuple with "
                             "length of 4")
        if not check_array_order(pre_filt, order="ascending"):
            raise ValueError("Input pre_filt must be in ascending order: %s"
                             % pre_filt)

    if remove_response_flag:
        # remove response
        if inventory is None:
            raise ValueError("Station information(inv) should be provided if"
                             "you want to remove instrument response")
        st.attach_response(inventory)
        if filter_flag:
            st.remove_response(output="DISP", pre_filt=pre_filt,
                               zero_mean=False, taper=False,
                               water_level=water_level)
        else:
            st.remove_response(output="DISP", zero_mean=False, taper=False)
    elif filter_flag:
        # Perform a frequency domain taper like during the response removal
        # just without an actual response...
        filter_stream(st, pre_filt)

    if filter_flag or remove_response_flag:
        # detrend, demean or taper
        st.detrend("linear")
        st.detrend("demean")
        st.taper(max_percentage=taper_percentage, type=taper_type)

    # resample
    if resample_flag:
        # interpolation
        if sampling_rate is None:
            raise ValueError("sampling rate should be provided if you set"
                             "resample_flag=True")

        if endtime is not None and starttime is not None:
            npts = int((endtime - starttime) * sampling_rate) + 1
            st = interpolate_stream(st, sampling_rate, starttime=starttime,
                                    npts=npts)
        else:
            # it doesn't matter if starttime is None or not, cause
            # obspy will handle this case
            st = interpolate_stream(st, sampling_rate, starttime=starttime)
    else:
        if starttime is not None and endtime is not None:
            # just cut
            st.trim(starttime, endtime)

    # rotate
    if rotate_flag:
        st = rotate_stream(st, event_latitude, event_longitude,
                           inventory=inventory, mode="ALL->RT",
                           sanity_check=sanity_check)

    # Convert to single precision to save space.
    for tr in st:
        tr.data = np.require(tr.data, dtype="float32")

    # transfer back to trace if input type is Trace
    if _is_trace:
        st = st[0]

    return st


if __name__ == '__main__':
    pass