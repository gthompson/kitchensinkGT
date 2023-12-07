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
import os, glob, sys

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.use('TkAgg')
sys.path.append('/Users/thompsong/src/volcanoObsPy/LIB')
from obspy import read, Stream, Trace
from obspy.core.utcdatetime import UTCDateTime
from converttime import datenum
import libseisGT as ls
import USF_instrument_responses as USF
import metrics
plt.close('all')
        
        
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
    ls.fix_seed_band_code(st)                    
    #print('all data:')
    #print(st.__str__(extended=True))
    return st

            
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

            USF.correctUSFstations(st)
            metrics.peak_amplitudes(st, eventdir)
            ls.plot_stream_types(st, eventdir)
        
            st.write(picklefile)
            
        #a = input('any key to continue')
    for tr in st:
         print(tr.stats)
         print(tr.stats.history)
         print(tr.stats.filter)
