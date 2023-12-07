# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#%reset
dropboxpath = '/Users/thompsong/Dropbox/PROFESSIONAL/RESEARCH/3 Project Documents/201602 Rocket Seismology/db/wfdata/rocketevents'
import os, glob
from obspy import read, Stream
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/thompsong/src/volcanoObsPy/LIB')
import libseisGT as ls
plt.close('all')
from obspy.core.utcdatetime import UTCDateTime as dt

def datenum(d):
    #print(d)
    #print([method for method in dir(d) if callable(getattr(d, method))])
    dnum = 366 + d.toordinal() + d.hour/24 + d.minute/3600 + d.second/86400 + d.microsecond/86400000000
    return dnum


def tellme(s):
    print(s)
    plt.title(s, fontsize=16)
    plt.draw()

os.chdir(dropboxpath)
eventdirs = sorted(glob.glob("20*"))
for eventdir in eventdirs:


    masterfile = os.path.join(eventdir, 'masterwaveformdata.mseed')
    if os.path.exists(masterfile):
        st = read(masterfile)
    else:

        # Load any seismo-acoustic data we can find in the event folder
        st = Stream()

        for root, dirs, files in os.walk(eventdir, topdown=False):
    
            for name in files:
                thisfullpath = os.path.join(root, name)
                if thisfullpath.lower().find('soh')==-1 and thisfullpath.lower().find('txt')==-1:
                    print(thisfullpath)
                    try:
                        this_st = read(thisfullpath)
                    except:
                        pass
                    else:
                        for tr in this_st:
                            if tr.stats.sampling_rate>=50:
                                tr.detrend()
                                if tr.stats.channel[1]=='H': # a seismic channel
                                    # trace is in counts. there are 0.3 counts/ (nm/s).
                                    tr.data = tr.data / 0.3 * 1e-9 # now in m/s
                                st.append(tr)
            
        # Master the data we have found for this event
        if st:
            for tr in st:
                if tr.stats.station == 'CARL0':
                    tr.stats.station = 'BCHH'
            #st.merge()
            #if st.select(station='BCHH') and st.select(station='CARL0'):
                #stgood = Stream()
                #for tr in st:
                    #if tr.stats.station != 'CARL0':
                        #stgood.append(tr)
                #st = stgood 
    
            # Save data
            print(st.__str__(extended=True))
            st.write(masterfile)      
            
    if st: 
        #plt.close('all')
        #ylabels = []
        #for tr in st:
        #    ylabels.append(tr.id)
        #ls.mulplt(st, 'UTC', ylabels)
        #plt.show()
        
        #a = input('any ley')
        plt.close('all')
        # select start and end points
        N = len(st)
        fig = plt.figure(figsize=[20,20])
        st.detrend()
        min_stime, max_stime, min_etime, max_etime = ls.Stream_min_starttime(st)
        st.trim(starttime=min_stime, endtime=max_etime, pad=True, fill_value=0)
        for c in range(N):
            tr = st[c]
            ax = fig.add_subplot(N, 1, N-c)
            ax.plot(tr.times("matplotlib"), tr.data, "b-")
            ax.xaxis_date()
            ax.set_ylabel(tr.id, rotation=0)
            ax.set_xlim(datenum(min_stime), datenum(max_etime))
        fig.autofmt_xdate()
        fig.show()
       
        #st.trim(starttime = x[0], endtime = x[1])

        #tellme('Now do a nested zoom, click to begin')
        #fig.waitforbuttonpress()

        got_two_points = False
        while not got_two_points:
            tellme('Select two corners of zoom, middle mouse button to finish')
            pts = fig.ginput(2, timeout=-1)
            if len(pts) < 2:
                print("user did not select two points")
            else:
                (x0, y0), (x1, y1) = pts
                xmin, xmax = sorted([x0, x1])
                ymin, ymax = sorted([y0, y1])
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
                got_two_points = True
            

        tellme('All Done!')
        plt.close('all')
        
        stZ = st.copy().select(channel="[ESBH]HZ")
        print('Synthetic 3-component vector amplitudes:')
        print('station, PGV (m/s), PGA (m/s-2)')
        for tr in stZ:
            thisID = tr.id[:-1]
            st3c = st.copy().select(id = '%s[ENZ12RT]' % thisID)
            if len(st3c)==3:
                mv = ls.max_3c(st3c)
                ma = ls.max_3c(st3c.differentiate())
                print('%s, %4.2e, %4.2e' % (thisID, mv[0], ma[0])) 
                tr.stats['PGV'] = mv[0]   
                tr.stats['PGA'] = ma[0]        
        
        # Plot velocity seismogram
        vpngfile = os.path.join(eventdir, 'velocity_seismogram.png')
        if not os.path.exists(vpngfile):
            st.plot(equal_scale=False, outfile=vpngfile)
         
    
        # Plot acceleration seismogram
        st.differentiate()
        apngfile = os.path.join(eventdir, 'acceleration_seismogram.png')
        if not os.path.exists(apngfile):
            st.plot(equal_scale=False, outfile=apngfile)

            
        a = input('any key to continue')
            
    # Note that this workflow does not correct the infrasound data channels: they are still in counts. 
    # And PGV and PGA values for infrasound channels are meaningless