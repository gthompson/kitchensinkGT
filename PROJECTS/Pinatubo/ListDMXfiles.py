#!/usr/bin/env python
# coding: utf-8

# Do all library imports first
import obspy.core as op
import glob
import matplotlib.pyplot as plt
import os
import shutil
import pandas as pd

# Set up all paths that will be used
os.chdir('D:\Dropbox\DATA\Pinatubo')
net = 'XB' # assigned by Gale Cox. 1R is code for KSC.
seisanDBname = 'PNTBO'
WAVEFORM_DIR = 'WAVEFORMS'

print('Paths setup okay')


# Loop over all files
alldirs = sorted(glob.glob(os.path.join(WAVEFORM_DIR, '*')))
lod = []
list_of_stations = []
for thisdir in alldirs:
    allDMXfiles = sorted(glob.glob(os.path.join(thisdir, '*.DMX')))
    for DMXfile in allDMXfiles:
        thisd = {'DMXfile':os.path.basename(DMXfile), 'time':0, 'nTraces':0, 'Fs':0, 'npts':0, 'WAV':''}
        try:
            st = op.read(DMXfile, headonly=True)
            for tr in st:
                list_of_stations.append(tr.stats.station)    
        except:
            print(os.path.basename(DMXfile), 0, 0, 0, 0)
        else:
            thisd['time']=st[0].stats.starttime
            thisd['nTraces']=len(st)
            thisd['Fs']=st[0].stats.sampling_rate
            thisd['npts']=st[0].stats.npts
            thisd['WAV']="%sM.%s_%03d" % (st[0].stats.starttime.strftime('%Y-%m-%d-%H%M-%S'), seisanDBname, len(st))
            #1995-01-23-1230-20M.BERGE_013
        print(thisd)
        lod.append(thisd)
df = pd.DataFrame(lod)
df.to_csv('ListDMXfiles_v2.csv')


import numpy as np
plt.figure()
ax = plt.subplot(111)
dtime=[]
for i,row in df.iterrows():
    if row['nTraces']>0:
        dtime.append(row['time'].datetime)
ax.plot(dtime, np.cumsum(np.ones(len(dtime))))
ax.xaxis_date()

plt.show()


print('number of DMX files = %d' % len(df.index))


list_of_stations = sorted(set(list_of_stations))
print('List of stations:\n',list_of_stations)


def fix_traceid(sta):
    if sta=='IRIG':
        return '.IRIG..'
    newchan = 'EH' + sta[-1]
    newsta = sta[:-1]
    loc = '00'
    
    traceid = '%s.%s.%s.%s' % (net, newsta, loc, newchan)
    return traceid
for sta in list_of_stations:
    print(sta, fix_traceid(sta))


#station_uptime = pd.DataFrame(columns=list_of_stations)
list_of_stations = ['IRIG', 'BUGZ', 'BURZ', 'CABN', 'CABZ', 'CRWZ', 'DONZ', 'FNGZ', 'GRNZ', 'PI2Z',
                    'PIEZ', 'PPOE', 'PPON', 'PPOZ', 'QADZ', 'UBOZ']
alldirs = sorted(glob.glob(os.path.join(WAVEFORM_DIR, '*')))
lod = []
for thisdir in alldirs:
    allDMXfiles = sorted(glob.glob(os.path.join(thisdir, '*.DMX')))
    for DMXfile in allDMXfiles:
        print(DMXfile)
        thisd = dict()
        for sta in list_of_stations:
            thisd[sta]=False
        try:
            st = op.read(DMXfile, headonly=True)
            for tr in st:
                thisd[tr.stats.station]=True  
        except:
            print(os.path.basename(DMXfile), 0, 0, 0, 0)
        #print(thisd)
        lod.append(thisd)
station_uptime = pd.DataFrame(lod)
station_uptime.to_csv('station_uptime.csv')
df = pd.concat([df.reset_index(drop=True),station_uptime.reset_index(drop=True)], axis=1)
df.to_csv('indexDMXfiles.csv')


print(station_uptime)


dfsum = station_uptime.sum(axis=0)
print(dfsum)


# Find out ondate and offdates for each station - INCOMPLETE
df = pd.read_csv('ListDMXfiles.csv')
station_uptime['time']=df['time']




