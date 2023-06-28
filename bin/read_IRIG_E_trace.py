#!/usr/bin/env python
# coding: utf-8
# based on https://safran-navigation-timing.com/manuals/SS/Content/NC_and_SS/Com/Topics/APPENDIX/IRIGe.htm
# looks like:
#     each element is 0.1 s or 10 samples
#     binary 0 is 20 ms, or 2 samples
#     binary 1 is 50 ms, or 5 samples
#     reference marker is 80 ms, or 8 samples, and should occur every 1-s or 100 samples
#     frame is 10 s. each sub-frame is 1 s.
#     subframes:
#               1: seconds
#               2: minutes
#               3: hours
#               4-5: days
#               6: sync status
#               7: year
#               9-10: binary seconds on this day
import os, sys
import obspy
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
import numpy as np
from numpy.random import rand
import datetime
plt.rcParams['figure.figsize'] = [16, 1]

verbose = True

def island_cumsum_vectorized(a):
    a_ext = np.concatenate(( [0], a, [0] ))
    idx = np.flatnonzero(a_ext[1:] != a_ext[:-1])
    a_ext[1:][idx[1::2]] = idx[::2] - idx[1::2]
    return a_ext.cumsum()[1:-1]

def detect_repeating_values_alt(a, n):
    b = island_cumsum_vectorized(a)
    if np.max(b)==n:
        return True
    else:
        return False

def detect_repeating_values(a, n):
    #print(a, n)
    if len(a)<n:
        return np.array([])
    N = n-1
    m = a[:-1]==a[1:]
    return np.flatnonzero(np.convolve(m,np.ones(N, dtype=int))==N)-N+1  # indices of where repeating digits found  

def detect_logical1_alt(a):
    return detect_repeating_values_alt(a, 5) 

def detect_logical1(a):
    d = detect_repeating_values(a, 5)
    if d.size==0:
        return 0
    else:
        return 1
    if d.any():
        return 1
    else:
        return 0
    return 

def detect_logical0(a):
    d = detect_repeating_values(a, 2)
    if d.any():
        return 1
    else:
        return 0
    return 

def detect_reference_markers(a):
    #print(a)
    # reference markers are 80 ms long, which means 8 samples at 100 Hz
    b = detect_repeating_values(a, 8)
    c = b.copy()
    for i in range(len(b)-1):
        if b[i+1]==b[i]+1:
            np.delete(c, i+1)
    return c
        

def detect_consecutive_reference_markers(a):
    # double reference marker marks start of frame. on time is the start of the second reference marker.
    # detect 8 high, 2 low, then 8 more high
    b = detect_reference_markers(a)
    print('Reference markers start at: ', b)
    c = []
    for i in range(1, len(b)):
        if b[i]==b[i-1]+10:
            c.append(b[i])
    return np.array(c)

def sample_rates(a):
    fs = []
    for i in range(1,len(a)):
        fs.append((a[i]-a[i-1])/10.0)
    return fs

def tracedata2logical(data):
    x = data/2048.0 # try to find 8 high samples in a row
    x = x + rand(1,len(x))/1e3
    x[x>0.5]=1.0    
    return x[0]

def get_bcd60(subseq):
    total = 0
    if subseq[0]==1: 
        total += 1
    if subseq[1]==1: 
        total += 2
    if subseq[2]==1: 
        total += 4 
    if subseq[3]==1: 
        total += 8
    if subseq[4]==1: # Dummy
        total += 0        
    if subseq[5]==1: 
        total += 10
    if subseq[6]==1:    
        total += 20
    if subseq[7]==1:
        total += 40
    if len(subseq)>8:
        if subseq[8]==1:
            total += 80       
    return total

def get_bcddays(seq):
    days = 0
    if seq[4][1]==1: 
        days += 200
    if seq[4][0]==1: 
        days += 100
    days += get_bcd60(seq[3])
    return days

def get_syncstatus(subseq):
    if subseq[4]==1: 
        return 1
    else:
        return 0
    
def get_binarysecs(seq):
    secs = 0
    for i in range(9):    
        if seq[8][i]==1: 
            secs += np.power(2,i)
        if seq[9][i]==1: 
            secs += np.power(2,i+9) 
    return secs
      
        
def read_frame(tr, si, ei):
    print('Frame from indices %d to %d' % (si, ei))
    s = tr.stats.starttime
    stime = s + si * tr.stats.delta
    print('Frame start time is ',stime)
    etime = s + ei * tr.stats.delta    
    #plt.plot(segment.times(), x, '.-')
    print('\nPlot of frame (10-s long):')
    tr.plot(starttime=stime, endtime=etime)
    
    x = tracedata2logical(tr.data[si:ei])
    refmarks = detect_reference_markers(x)
    print(refmarks)
    plt.plot(list(range(len(x))), x, '-')
    for refmark in refmarks:
        plt.plot(refmark, 0, 'ro')
    plt.show()
            
    # build sequence
    sequence = []
    subframeNames = ['seconds', 'minutes', 'hours', 'days1', 'days2', 'sync_status', 'years', 'dummy', 
                     'binsecs1', 'binsecs2']
    for ii, refmark in enumerate(refmarks[0:-1]):
        subsequence = []
        if verbose:
            print('P%d: Sub-Frame %s from indices %d to %d' % (ii, subframeNames[ii], refmark, refmarks[ii+1]))
            #plt.plot(x[refmark+0:refmark+90], '.-')
            #plt.show()
        r = np.arange(10,refmarks[ii+1] - refmark, 10)
        print('r=',r)
        for start_i in r:
            xis = refmark + start_i
            xie = xis + 10
            
            print('length(x)=',len(x))
            print('subsetting x from indices %d to %d' % (xis,xie))
            xsubset = x[xis:xie]
            print('xsubset=',xsubset)
            res = detect_logical1(xsubset)
            print('result=',res)
            subsequence.append(res)
            #if verbose:
            #    print(np.abs(np.round(xsubset)), res)
        sequence.append(subsequence)
    print(sequence)
    secs = get_bcd60(sequence[0])
    mins = get_bcd60(sequence[1])
    hrs = get_bcd60(sequence[2])
    days = get_bcddays(sequence)
    years = get_bcd60(sequence[6])
    if years==0:
        years = tr.stats.starttime.year
    print('synced ?: ',get_syncstatus(sequence[5]))
    print('Y,Mo,D,H,Mi,S: ',years,days,hrs,mins,secs)
    binarysecs = get_binarysecs(sequence)
    print('binarysecs = ',binarysecs)
    print('ObsPy Frame start time is ',stime) 
    irigdatetime = datetime.datetime(years, 1, 1) + datetime.timedelta(days - 1 + hrs/24 + mins/1440 + secs/86400)
    print('IRIG frame start time is',irigdatetime)
    print('\n')
    return 0

# Main program
print(sys.argv)
if len(sys.argv)>1:
    seismicfile=sys.argv[1]
    if not os.path.isfile(seismicfile):
        raise ValueError('%s not found.' % seismicfile)
else:
    root = tk.Tk()
    root.withdraw()
    print('Select file with open file dialog')
    seismicfile = filedialog.askopenfilename()
    print('You chose: %s' % seismicfile)
try:
    st = obspy.read(seismicfile)
except:
    raise ValueError('obspy cannot read %s' % seismicfile)
for ti, tr in enumerate(st):
    print('Trace %d of %d: ID=%s' % (ti+1, len(st), tr.id))
    if tr.stats.station=='IRIG':
        tr.plot()
        print(tr.data)
        if 'DMX' in seismicfile:
            tr.data = tr.data - 2048.0
        a = tr.data.copy()
        x = tracedata2logical(a)
        frame_start_indices = detect_consecutive_reference_markers(x)
        print('Frames start at: ', frame_start_indices)
        fs = sample_rates(frame_start_indices)
        print('Sample rates (Hz) of each frame: ', fs)
        for index in range(len(frame_start_indices)-1):
            si = frame_start_indices[index]
            ei = frame_start_indices[index+1] #-1
            print('*************************')
            print('Examining COMPLETE FRAME %d of %d' % (index+1, len(frame_start_indices)-1))
            rtncode = read_frame(tr, si, ei )
            if index < len(frame_start_indices)-2: # Hit ENTER to see next frame
                anykey = input('<ENTER> to continue')


# In[ ]:




