# -*- coding: utf-8 -*-
"""
Process each rocket event.

Things to add:
    1) Wrap a website around these plots and PG values. Add spectrograms. Make it easy to navigate from one event to next.
    2) Run a detector like Miami Lakes? Put events into Seisan DB.
    3) Use moving window x-correlations in arrays to determine best azimuth and apparent speed and plot them as a function of time
    4) Use above to guess when rocket takes off. This could help Doppler-based methods.
    5) ASL. Try on seismic stations only, and infrasound stations only.
    6) Examine seismic wave types -> GISMO?
    7) Test that equation that Cassandra used for Pavlof infrasound to seismic

    * Check infrasound calibrations/instruments used
    * extract useful functions into libseisGT
"""
#%reset
dropboxpath = '/Users/thompsong/Dropbox/PROFESSIONAL/RESEARCH/3 Project Documents/201602 Rocket Seismology/db/wfdata/rocketevents'
import os, glob
from obspy import read, Stream, Trace
import numpy as np
#import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/thompsong/src/volcanoObsPy/LIB')
import libseisGT as ls
plt.close('all')
from obspy.core.utcdatetime import UTCDateTime



def datenum(d):
    #print(d)
    #print([method for method in dir(d) if callable(getattr(d, method))])
    dnum = 366 + d.toordinal() + d.hour/24 + d.minute/3600 + d.second/86400 + d.microsecond/86400000000
    return dnum



def fix_seed_band_code(st):
    for tr in st:
        sr = tr.stats.sampling_rate
        if sr > 1 and sr < 10:
            bc = 'M'
        if sr >= 10 and sr < 80:
            bc = 'B'
        if sr >= 80:
            bc = 'H'
        if not bc == tr.stats.channel[0]:
            tr.stats.channel = bc + tr.stats.channel[1:]
            add_to_trace_history(tr, 'bandcode_fixed') 
        
        
def load_all_data(eventdir):
    # Load any seismo-acoustic data we can find in the event folder
    st = Stream()

    for root, dirs, files in os.walk(eventdir, topdown=False):
        if 'obsolete' in root.lower():
            continue  
            
        for name in files:
            if '.DS_store' in name:
                continue
            thisfullpath = os.path.join(root, name)
            print(root, name)
            if thisfullpath.lower().find('soh')==-1 and thisfullpath.lower().find('txt')==-1 and thisfullpath.lower().find('png')==-1:
                print('Attempting to read %s ' % thisfullpath)
                try:
                    this_st = read(thisfullpath)
                    print('- succeeded')
                except:
                    print('- failed')
                else:
                    for tr in this_st:
                        if tr.stats.station == 'CARL0':
                            tr.stats.station = 'BCHH' 
                            add_to_trace_history(tr, 'station_fixed') 
                        if not tr.stats.channel[0]=='L' and tr.stats.sampling_rate > 10:                          
                            st.append(tr)
    fix_seed_band_code(st)                    
    #print('all data:')
    #print(st.__str__(extended=True))
    return st


def add_to_trace_history(tr, str):
    if not 'history' in tr.stats:
        tr.stats['history'] = list()
    if not str in tr.stats.history:
        tr.stats.history.append(str)

def update_trace_filter(tr, filtertype, freq, causal):
    if not filter in tr.stats:
        tr.stats['filter'] = {'freqmin':0, 'freqmax':tr.stats.sampling_rate/2, 'causal': False}
    if filtertype == 'highpass':    
        tr.stats.filter["freqmin"] = max([freq, tr.stats.filter["freqmin"]])
    if filtertype == 'bandpass':
        tr.stats.filter["freqmin"] = max([freq[0], tr.stats.filter["freqmin"]]) 
        tr.stats.filter["freqmax"] = min([freq[1], tr.stats.filter["freqmax"]])
    if filtertype == 'lowpass':
        tr.stats.filter["freqmax"] = min([freq, tr.stats.filter["freqmax"]])
    tr.stats.filter['causal'] = causal

            
def input_trim_times(st):        
    min_stime, max_stime, min_etime, max_etime = ls.Stream_min_starttime(st)
    sy = min_stime.year
    sm = min_stime.month
    sd = min_stime.day
    if max_etime.day != min_stime.day:
        symd = input('Start year/month/day (default: %d/%d/%d): ' % (sy,sm,sd))
        if symd:
            sy, sm, sd = symd.split('/')
            sy = int(sy)
            sm = int(sm)
            sd = int(sd)
    shour = int(input('Start hour: '))
    smin = int(input('Start minute: '))

    ey = sy
    em = sm
    ed = sd                
    if max_etime.day != min_stime.day:
        eymd = input('End year/month/day: (default: %d/%d/%d): ' % (ey,em,ed))
        if eymd:
            ey, em, ed = eymd.split('/')
            ey = int(ey)
            em = int(em)
            ed = int(ed)           
    ehour = int(input('End hour: '))
    emin = int(input('End minute: ') )

    manstime = UTCDateTime(sy, sm, sd, shour, smin, 0)
    manetime = UTCDateTime(ey, em, ed, ehour, emin, 0)
        
    return (manstime, manetime)     


def manual_zoom(st, eventdir):            
    still_zooming = True
    while still_zooming:
        st.plot(equal_scale=False, outfile=os.path.join(eventdir, 'zooming.png'))
        min_stime, max_stime, min_etime, max_etime = ls.Stream_min_starttime(st)
        print('Stream start: %s' % min_stime)
        print('Stream end:   %s' % max_etime)
        yn = input('Do you want manually trim times further? ')
        if yn.lower()[0]=='y':
            newstime, newetime = input_trim_times(st)
            st.trim(starttime=newstime, endtime=newetime, pad=False, fill_value=0)
        else:
            still_zooming = False



def smart_merge_traces(trace_pair):
    """
    Clever way to merge overlapping traces. Uses all non-zero data values from both.
    """
    this_tr = trace_pair[0] 
    other_tr = trace_pair[1]

    error_flag = False


    if not (this_tr.id == other_tr.id):
        print('Different trace IDs. Cannot merge.')
        error_flag = True

    if not (this_tr.stats.sampling_rate == other_tr.stats.sampling_rate):
        print('Different sampling rates. Cannot merge.')
        error_flag = True

    if (abs(this_tr.stats.starttime - other_tr.stats.starttime) > this_tr.stats.delta/4):
        print('Different start times. Cannot merge.')
        error_flag = True

    if (abs(this_tr.stats.endtime - other_tr.stats.endtime) > this_tr.stats.delta/4):
        print('Different end times. Cannot merge.')
        error_flag = True

    if error_flag: # traces incompatible, so return the trace with the most non-zero values
        this_good = np.count_nonzero(this_tr.data)
        #print(this_tr.stats)
        other_good = np.count_nonzero(other_tr.data)
        #print(other_tr.stats)
        if other_good > this_good:
            return other_tr
        else:
            return this_tr

    else: # things are good
        indices = np.where(other_tr.data == 0)
        other_tr.data[indices] = this_tr.data[indices]
        return other_tr

            
def smart_merge(st):
    # need to loop over st and find traces with same ids
    ##### GOT HERE
    newst = Stream()
    all_ids = []
    for tr in st:
        if not tr.id in all_ids:
            all_ids.append(tr.id)

    for this_id in all_ids: # loop over all nsl combinations
        these_traces = st.copy().select(id=this_id).sort() # find all traces for this nsl combination
        
        # remove duplicates
        traces_to_remove = []
        for c in range(len(these_traces)-1):
            s0 = these_traces[c].stats
            s1 = these_traces[c].stats
            if s0.starttime == s1.starttime and s0.endtime == s1.endtime and s0.sampling_rate == s1.sampling_rate:
                traces_to_remove.append(c)
                
        print(these_traces)
        print(traces_to_remove)
        if traces_to_remove:
            for c in traces_to_remove:
                these_traces.remove(these_traces[c])
                
        if len(these_traces)==1: # if only 1 trace, append it, and go to next trace id
            newst.append(these_traces[0]) 
            continue
        
        # must have more than 1 trace
        try: # try regular merge now duplicates removed
            merged_trace = these_traces.copy().merge()
            print('- regular merge of these traces success')
        except:
            print('- regular merge of these traces failed')   
            # need to try merging traces in pairs instead
            N = len(these_traces)
            these_traces.sort() # sort the traces
            for c in range(N-1): # loop over traces in pairs
                appended = False
                
                # choose pair
                if c==0:
                    trace_pair = these_traces[0:2]
                else:
                    trace_pair = Stream(traces=[merged_trace, these_traces[c+1] ] )
                                        
                # merge these two traces together    
                try: # standard merge
                    merged_trace = trace_pair.copy().merge()
                    print('- regular merge of trace pair success')
                except: # smart merge
                    print('- regular merge of trace pair failed')
                    try:
                        min_stime, max_stime, min_etime, max_etime = ls.Stream_min_starttime(trace_pair)
                        trace_pair.trim(starttime=min_stime, endtime=max_etime, pad=True, fill_value=0)
                        merged_trace = Stream.append(smart_merge_traces(trace_pair)) # this is a trace, not a Stream
                    except:
                        print('- smart_merge of trace pair failed')
                        
            # we have looped over all pairs and merged_trace should now contain everything
            # we should only have 1 trace in merged_trace
            print(merged_trace)
            if len(merged_trace)==1:
                try:
                    newst.append(merged_trace[0])
                    appended = True
                except:
                    pass
                
            if not appended:
                print('\n\nTrace conflict\n')
                trace_pair.plot()
                for c in range(len(trace_pair)):
                    print(c, trace_pair[c])
                choice = int(input('Keep which trace ? '))
                newst.append(trace_pair[choice])  
                appended = True                
                             
    return newst






        
def centaur(inputVoltageRange):
    countsPerVolt = 0.4e6 * 40/inputVoltageRange;
    return countsPerVolt

def trillium():
    voltsPerMS = 750; # V / (m/s)
    return voltsPerMS

def infraBSU(HgInThisSensor=0.5):
    # model 0.5" is default
    # 0.1 - 40 Hz flat
    oneInchHg2Pa = 3386.4;
    linearRangeInPa = oneInchHg2Pa * HgInThisSensor;
    selfNoisePa = 5.47e-3;
    voltsPerPa = 46e-6; # from infraBSU quick start guide
    return voltsPerPa

def ChaparralM25():
    # 36 V p2p
    # 0.1 - 200 Hz flat
    selfNoisePa = 3e-3;
    voltsPerPaHighGain = 2.0; # 18 Pa full scale. 
    voltsPerPaLowGain = 0.4; # 90 Pa full scale. recommended for 24-bit digitizers. 
    voltsPerPaVMod = 0.4 * 90/720; # 720 Pa full scale.
    # Volcano mod reduces sensitivity further.
    return voltsPerPaLowGain

countsPerMS = centaur(40.0) * trillium()
countsPerPa40 = centaur(40.0) * infraBSU(0.5)
countsPerPa1 = centaur(1.0) * infraBSU(0.5)
countsPerPaChap = centaur(40.0) * ChaparralM25()

def correctUSFstations(st):
    for tr in st:
        if tr.stats.sampling_rate>=50:
            #tr.detrend()
            if tr.stats.channel[1]=='H': # a seismic velocity high-gain channel. L for low gain, N for accelerometer
                # trace is in counts. there are 0.3 counts/ (nm/s).
                #tr.data = tr.data / 0.3 * 1e-9 # now in m/s
                calib = countsPerMS
                units = 'm/s'
                                
            if tr.stats.channel[1]=='D': # infraBSU channel?
                #trtr.data = tr.data / countsPerPa40
                # Assumes 0.5" infraBSU sensors at 40V p2p FS
                # But could be 1" or 5", could be 1V p2p FS or could be Chaparral M-25
                units = 'Pa'
                if tr.stats.station=='BCHH1' and tr.stats.channel[2]=='1': # Chaparral M-25
                    calib = countsPerPaChap
                elif tr.stats.station=='BCHH' or tr.stats.station=='BCHH1' or tr.stats.station[0:3]=='SAK' or tr.stats.network=='NU':
                    calib = countsPerPa40
                else:
                    calib = countsPerPa1
        
            tr.data = tr.data / calib
            tr.stats['calib'] = calib
            tr.stats['units'] = units
            add_to_trace_history(tr, 'calibration_applied')        
        



def clean_stream(st, taper_fraction=0.05, freqmin=0.05, causal=True):
    """
    Clean Stream object in place.
    clean_stream(st, taper_fraction=0.05, freqmin=0.05, causal=True)
    detrend, taper and high pass
    """
    if causal:
        causalstr = 'causal'
    else:
        causalstr = 'acausal'
    for tr in st:
        if not 'history' in tr.stats:
            tr.stats['history'] = list()
        try:
            if not 'detrended' in tr.stats.history:
                tr.detrend()
                add_to_trace_history(tr, 'detrended')
        except:
            pass
        try:
            if not 'tapered' in tr.stats.history:
                tr.taper(max_percentage=taper_fraction, type="hann") 
                add_to_trace_history(tr, 'tapered')
        except:
            pass
        try:
            filtertype = 'highpass'
            tr.filter(filtertype, freq=0.05, corners=2, zerophase=causal)
            update_trace_filter(tr, filtertype, freqmin, causal)
            add_to_trace_history(tr, filtertype)
        except:
            pass


def peak_amplitudes(st, eventdir):   
    eventbase = eventdir.split('/')[-1]
    
    clean_stream(st, taper_fraction=0.05, freqmin=0.05, causal=True)
               
    # velocity, displacement, acceleration seismograms
    stV = st.copy().select(channel='[ESBH]H?') 
    stD = stV.copy().integrate()
    for tr in stD:
        add_to_trace_history(tr, 'integrated')    
    stA = stV.copy().differentiate()
    for tr in stD:
        add_to_trace_history(tr, 'differentiated') 
     
    # Seismic vector data
    peakseismofile = os.path.join(eventdir, 'summary_seismic_3c.csv')
    fptr = open(peakseismofile, 'w')
    fptr.write('eventdir, id, PGD, PGV, PGA, calib, units\n')    
    stZ = stV.copy().select(channel="[ESBH]HZ")
    for tr in stZ:
        thisID = tr.id[:-1]
        st3cV = stV.copy().select(id = '%s[ENZ12RT]' % thisID)
        if len(st3cV)==3:
           st3cD = stD.copy().select(id = '%s[ENZ12RT]' % thisID)
           st3cA = stA.copy().select(id = '%s[ENZ12RT]' % thisID)
           md = ls.max_3c(st3cD)
           mv = ls.max_3c(st3cV)
           ma = ls.max_3c(st3cA)
           fptr.write('%s, %s, %4.2e, %4.2e, %4.2e, %4.2e, %s\n' % (eventbase, thisID, md[0], mv[0], ma[0], tr.stats.calib, tr.stats.units)) 
           tr.stats['PGD'] = md[0]
           tr.stats['PGV'] = mv[0]   
           tr.stats['PGA'] = ma[0]
    fptr.close()
    
    # Seismic 1-c data
    peakseismo1cfile = os.path.join(eventdir, 'summary_seismic_1c.csv')
    fptr = open(peakseismo1cfile, 'w')
    fptr.write('eventdir, id, PGD, PGV, PGA, calib, units\n')
    for c in range(len(stV)):
        md = max(abs(stD[c].data))
        mv = max(abs(stV[c].data))
        ma = max(abs(stA[c].data)) 
        fptr.write('%s, %s, %4.2e, %4.2e, %4.2e, %4.2e, %s\n' % (eventbase, stV[c].id, md, mv, ma, stV[c].stats.calib, stV[c].stats.units)) 
    fptr.close()    
            
    # Infrasound data
    peakinfrafile = os.path.join(eventdir, 'summary_infrasound.csv')
    fptr = open(peakinfrafile, 'w')
    fptr.write('eventdir, id, broadband, bandpassed, calib, units\n')
    stP = st.copy().select(channel="[ESBH]D?")
    stPF = stP.copy().filter('bandpass', freqmin=1.0, freqmax=20.0, corners=2, zerophase=True)    
    for c in range(len(stP)):
        mb = max(abs(stP[c].data))
        mn = max(abs(stPF[c].data)) 
        fptr.write('%s, %s, %4.2e, %4.2e, %4.2e, %s\n' % (eventbase, stP[c].id, mb, mn, stP[c].stats.calib, stP[c].stats.units)) 
    fptr.close()
    
    
    # Plot displacement seismogram
    if stD:
        dpngfile = os.path.join(eventdir, 'seismogram_D.png')
        stD.plot(equal_scale=False, outfile=dpngfile)    
    
    # Plot velocity seismogram
    if stV:
        vpngfile = os.path.join(eventdir, 'seismogram_V.png')
        stV.plot(equal_scale=False, outfile=vpngfile)
         
    # Plot acceleration seismogram
    if stA:
        apngfile = os.path.join(eventdir, 'seismogram_A.png')
        stA.plot(equal_scale=False, outfile=apngfile)
        
    # Plot pressure acoustograms
    if stP:
        ppngfile = os.path.join(eventdir, 'seismogram_P.png')
        stP.plot(equal_scale=False, outfile=ppngfile)        
                     
                   

#########################
#### MAIN PROGRAM #######
#########################
os.chdir(dropboxpath)
eventdirs = sorted(glob.glob("20*"))
for eventdir in eventdirs:
    print('-' * 50)
    print(eventdir)
    
    if 'sonic' in eventdir:
        continue

    picklefile = os.path.join(eventdir, 'trimmed.pickle')
    trimmedfile = os.path.join(eventdir, 'trimmed.mseed')
    allfile = os.path.join(eventdir, 'all.mseed')
    if os.path.exists(picklefile):
        st = read(picklefile)    
    else:
        if os.path.exists(trimmedfile):
            st = read(trimmedfile)
        else:
            if os.path.exists(allfile):
                st = read(allfile)
            else:
                st = load_all_data(eventdir)
            if st:
                st.plot(equal_scale=False, outfile=os.path.join(eventdir, 'all.png'))
                st.write(allfile)
                manual_zoom(st, eventdir)
                st.write(trimmedfile)
                                
        if st:        
            try:
                st.merge()
                print('standard merge success')
            except:
                print('standard merge failed')
                st = smart_merge(st)

            correctUSFstations(st)
            peak_amplitudes(st, eventdir)
        
            st.write(picklefile)
            
        #a = input('any key to continue')
    for tr in st:
         print(tr.stats)
         print(tr.stats.history)
         print(tr.stats.filter)
