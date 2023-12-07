#!/usr/bin/env python
# coding: utf-8

# Glenn Thompson, Fall 2020 - May 12, 2021

# imports
import pandas as pd
import obspy.core
#import numpy as np
#from sys import exit
import os, sys
import glob
LIBpath = os.path.join( os.getenv('HOME'),'src','kitchensinkGT', 'LIB')
sys.path.append(LIBpath)
from libseisGT import inventory2traceid, get_FDSN_inventory, get_FDSN_Stream, removeInstrumentResponse, detectEvent
from libgeoGT import km2deg
from metrics import eventStatistics

### SCAFFOLD. To do: replace eventStatistics with other routines in metrics.py that measure more metrics.

# constants
TOPDIR =  os.path.join(os.getenv('HOME'), 'Miami_Lakes2')
PRE_TRIGGER_SECS = 300
POST_TRIGGER_SECS = 300
LATITUDE = 25.911667
LONGITUDE = -80.325
NETWORK = "AM"
searchRadiusKm = 20
searchRadiusDeg = km2deg(searchRadiusKm)
TEST_MODE = True

from obspy.clients.fdsn import Client
fdsnClient = Client("RASPISHAKE")

# Read Steve's Excel file (after conversion to CSV) into a dataFrame, df
#excelfile = os.path.join(TOPDIR, 'MiamiBlastsMasterList.xlsx')
excelfile = 'MiamiBlastsMasterList.xlsx'
df = pd.read_excel(excelfile, sheet_name='All')
df.dropna(subset=['UTCdatetime'], inplace=True)
eventTime = list()
for index, row in df.iterrows():
    thisEventTime = obspy.core.utcdatetime.UTCDateTime(row['UTCdatetime'])
    eventTime.append(thisEventTime)
df['time'] = pd.Series(eventTime, index=df.index)
#df = df[["event_num", "time", "PPV_max", "Pa", "type", "AW-GW", "distance", "reduced_pressure", "comment"]]
df = df[["event_num", "time", "PPV_max", "Pa", "type" ]]
print(df)

if TEST_MODE:
    df = df.tail()

# Save this improved event CSV file
csvmasterfile = os.path.join(TOPDIR, 'MiamiBlastsMasterListPandas.csv')
df.to_csv(csvmasterfile, index=False)
# # Read the CSV file that has the time column added
# df = pd.read_csv(csvmasterfile)      

''' ------------ Main loop over events starts here --------------'''
for index, row in df.iterrows():
    yyyymmddhhmm = row['time'].strftime('%Y%m%d%H%M')
    stationXmlFile = os.path.join(TOPDIR, 'events', 'inventories', 'Inventory.%s.xml' % yyyymmddhhmm)
    inv = get_FDSN_inventory(fdsnClient, row['time'], stationXmlFile, NETWORK, LATITUDE, LONGITUDE, searchRadiusDeg, PRE_TRIGGER_SECS, POST_TRIGGER_SECS )
    if not inv:
        print('No inventory returned for event at ',row['time'])
        continue
            
    #inv.plot(projection='local', resolution='l')
    eventTime = list()

    # for sta in inv[0].stations:
    #     print(sta)
    
    # get list of unique trace ids
    trace_ids = inventory2traceid(inv)    

    # file paths to save
    rawfile = os.path.join(TOPDIR, 'events', 'mseed', 'raw', 'Event.%s.RAW.mseed' % yyyymmddhhmm)
    correctedfile = os.path.join(TOPDIR, 'events', 'mseed', 'corrected', 'Event.%s.CORRECTED.mseed' % yyyymmddhhmm)
    #pngfile = os.path.join(TOPDIR, 'events', 'plots', 'corrected', 'Event.%s.png' % yyyymmddhhmm)
    
    # load raw data and save to rawfile
    startt = row['time'] - PRE_TRIGGER_SECS
    endt = row['time'] + POST_TRIGGER_SECS    
    stR = get_FDSN_Stream(fdsnClient, trace_ids, rawfile, startt, endt )
    
    # detect events inside this Stream
    #stD = detect_network_event(st)
    # detect_network_event(st, minchans=None, threshon=3.5, threshoff=1.0, sta=0.5, lta=5.0, pad=0.0)
    
    # Load corrected data, or correct from raw data
    if os.path.exists(correctedfile):
        # Load corrected data from file
        stC = obspy.core.read(correctedfile)
    else:
        # Correct raw data & save to file
        stC = removeInstrumentResponse(stR, preFilter = (0.5, 0.75, 30.0, 45.0), outputType = "VEL")
"""        
    # IMPORTANT #
    # Now trim stC to boundaries of stD. We do this because correcting stD would result in artefacts within the measurement window. Better to use the longer time series for instrument correction, then trim and measure    
    stC.trim(starttime = stD[0].stats.starttime, endtime = stD[0].stats.endtime)

    # Seisan event name for this detected Stream
    isofmt = stC[0].stats.starttime.isoformat()
    seisanmseed = isofmt[0:10] + "-" + isofmt[11:13] + isofmt[14:16] + "-" + isofmt[17:19] + 'S.MILAK__%03d' % (len(st)+1)
    
    # Write detected event
    print('- writing %s' % seisanmseed)
    stC.write(os.path.join(TOPDIR, 'events', 'mseed', 'detected', seisanmseed), format='MSEED' )    
    
    # Plot corrected event
    pngfile = os.path.join(TOPDIR, 'events', 'plots', 'detected', seisanmseed + '.png')
    stC.plot(outfile = pngfile, equal_scale = False)
    
"""
    # Make some measurements & save to CSV file
    statsfile = os.path.join(TOPDIR, 'events', 'statistics', 'Event.%s.stats.csv' % yyyymmddhhmm)
    if not os.path.exists(statsfile):
        eventDF = eventStatistics(stC)
        eventDF.to_csv(statsfile, index=False)

        
# Update AllEventStats.csv 

STATSDIR = os.path.join(TOPDIR, 'events', 'statistics')
allDF = pd.DataFrame()
allEventFiles = sorted(glob.glob(os.path.join(STATSDIR, 'Event.2*.stats.csv') ))
numEvents = len(allEventFiles)
count = 0
for thisEventFile in allEventFiles:
    count = count + 1
    print("Processing %s (%d of %d)" % (thisEventFile, count, numEvents) )
    thisDF = pd.read_csv(thisEventFile)
    if not allDF.empty:
        allDF = allDF.append(thisDF)
    else:
        allDF = thisDF

allDF.to_csv(os.path.join(STATSDIR, 'AllEvents.stats.csv' ))

