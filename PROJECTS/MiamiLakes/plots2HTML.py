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
def plotAllStats(DF, selectedTraceID):

    net, sta, loc, chan = selectedTraceID.split('.')
    dates = list()
    for thistime in DF['time']:
        dates.append(UTCDateTime(thistime).datetime)

    fig = plt.figure(figsize=(20,10))
    ax = plt.subplot(121)
    ax.plot_date(dates, DF['peakamp'], '.')
    ax.set_title(selectedTraceID)
    ax.set_xlabel('Date')
    if chan[1] == 'H': # seismic
        ax.set_ylabel('Peak Amplitude (m/s)')
    else:
        ax.set_ylabel('Peak Amplitude (Pa)')
    for label in ax.get_xticklabels():
        label.set_rotation(30)
        label.set_horizontalalignment('right')

    ax = plt.subplot(122)
    ax.plot_date(dates, DF['energy'], '.')
    ax.set_title(selectedTraceID)
    ax.set_xlabel('Date')
    ax.set_ylabel('Relative Energy')
    for label in ax.get_xticklabels():
        label.set_rotation(30)
        label.set_horizontalalignment('right')

    plotPath = os.path.join('Figures', 'AllEvents_%s.png' % selectedTraceID) 
    fig.savefig(plotPath)
    plt.close()

    return(plotPath)

f = open('index.html','w')

message = """<html>\n"
<head>
</head>
<body>"""

uniqueTraceIDs = allDF['id'].unique()
count = 0
message += "<table border=1>\n"
message += "<tr><th>Index</th><th>traceID</th><th>numberOfEvents</th><th>thumbnail plots</th></tr>\n"
for thisID in uniqueTraceIDs:
    count += 1
    thisDF = allDF[allDF['id'] == thisID]
    print(count, thisID, len(thisDF.index) )
    plotPath = plotAllStats(thisDF, thisID)
    plotLink = "<a href='%s'> <img src='%s' height=200 /></a>" % (plotPath, plotPath)

    message += "<tr><td>%d</td><td>%s</td><td>%d</td><td>%s</td></tr>\n" % (count, thisID, len(thisDF.index), plotLink )

message += "</table>\n"
message += """</body>
</html>"""
f.write(message)
f.close()
