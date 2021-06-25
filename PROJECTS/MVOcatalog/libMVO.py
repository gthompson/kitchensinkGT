#!/usr/bin/env python
"""
function w = clean(w)
% CLEAN Clean up waveform object(s)
%   w = clean(w) will detrend waveforms (in a smart way, aware of NaN
%   values marking missing values), then use fillgaps to mark bounded NaNs
%   with linear-interpolated values, and also mark bounded NaNs with zeroes
%   (now okay, since trend has been removed, so no weird startup effects).
%   It then removes non-linear trends using a 20-s (0.05 Hz) highpass
%   filter.

% Glenn Thompson March 28, 2017

    for c=1:numel(w)
        
        if ~isempty(w(c))
            
            % remove spikes of length 1
            w(c) = medfilt1(w(c), 3); 

            % smart detrend
            w(c) = detrend(w(c));

            % % fill gaps to get rid of NaNs marking missing values, so we can filter
            w(c) = fillgaps(w(c), 'interp');

            % highpass data at 20s - to remove non-linear trends
            f = filterobject('h',0.05,2);
            w(c) = filtfilt(f,w(c));
        
        end
    end

end
"""

def check0andMinus1(liste):
    liste=list(liste)
    listStr=''.join(str(i) for i in liste)
    if  "000000000000" in listStr or "-1-1-1-1-1-1-1-1" in listStr :
        return False
    else:
        return True

def fix_trace_ids(tr): # add stuff here to correct sampling rate too
    # fix the network, channel and location
    network = 'MV'
    tr.stats['network']=network
    sta = tr.stats['station'].strip()
    chan = tr.stats['channel'].strip()
    if chan=='PRS' or chan=='APS':
        chan='BDF'
    else:
        if chan[0]=='A':
            if tr.stats['location'] == 'J':
                bandcode = 'S'
            #else:
                #bandcode = 'B'
        else:
            if chan[1]=='B':
                bandcode = 'B'
            else:
                bandcode = chan[0]
            instrumentcode = 'H'
            if len(chan)==2:
                orientationcode = tr.stats['location']
            else:
                orientationcode = chan[2]
            chan = bandcode + instrumentcode + orientationcode

    if chan[0]=='A':
        print(tr.stats)
        print(chan)
        sys.exit('bugger!')
    tr.stats['channel'] = chan
    tr.stats['location']='--'
    return tr

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
        

        



            
    
def reconstituteTrace(tr, taperFraction=0.05, filterType="bandpass", freq=[0.1, 20.0], corners=2, zerophase=True, inv=None):

    traceIsClipped = detectClipping(tr) # can add another function to interpolate clipped values
    
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
