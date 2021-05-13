#!/usr/bin/env python
# coding: utf-8

# Glenn Thompson, Feb 2021

# imports
import pandas as pd
import os
#import glob
import matplotlib.pyplot as plt
#import datetime
from obspy.core import UTCDateTime

STATSDIR = "eventStats"
allDF = pd.read_csv(os.path.join(STATSDIR, 'AllEvents.stats.csv' ))

uniqueTraceIDs = allDF['id'].unique()
count = 0
print('Index, traceID, numberOfEvents' )
for thisID in uniqueTraceIDs:
    count += 1
    thisDF = allDF[allDF['id'] == thisID]
    print(count, thisID, len(thisDF.index) )

tracenum = input('Enter integer corresponding to the traceID you want to plot: ? ')
selectedTraceID = uniqueTraceIDs[int(tracenum)-1]
print('OK, extracting all events for %s' % selectedTraceID)
#print(allDF)

subsetDF = allDF[allDF['id'] == selectedTraceID]
#dates = datetime.datetime.strptime(subsetDF['time'], '%Y%j')
#dates = UTCDateTime(subsetDF['time']).datetime
dates = list()
for thistime in subsetDF['time']:
    dates.append(UTCDateTime(thistime).datetime)
print('Number of corresponding events = %d' % len(subsetDF.index) )
#print(subsetDF)
net,sta,loc,chan = selectedTraceID.split('.')

fig = plt.figure()
ax = plt.subplot(111)
ax.plot_date(dates, subsetDF['peakamp'], 'o')
ax.set_title(selectedTraceID)
ax.set_xlabel('Date')
if chan[1] == 'H': # seismic
    ax.set_ylabel('Peak Amplitude (m/s)')
else:
    ax.set_ylabel('Peak Amplitude (Pa)')
for label in ax.get_xticklabels():
    label.set_rotation(30)
    label.set_horizontalalignment('right')
fig.savefig(os.path.join('Figures', '%s_peakamp.png' % selectedTraceID) )

fig = plt.figure()
ax = plt.subplot(111)
ax.plot_date(dates, subsetDF['energy'], 'o')
ax.set_title(selectedTraceID)
ax.set_xlabel('Date')
ax.set_ylabel('Relative Energy')
for label in ax.get_xticklabels():
    label.set_rotation(30)
    label.set_horizontalalignment('right')
fig.savefig(os.path.join('Figures', '%s_energy.png' % selectedTraceID) )


