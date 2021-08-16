#!/usr/bin/env python
# coding: utf-8

# Glenn Thompson, Fall 2020 - May 12, 2021

# imports
import pandas as pd
import obspy.core
#import numpy as np
#from sys import exit
import os
import glob
import sys
LIBpath = os.path.join( os.getenv('HOME'),'src','kitchensinkGT', 'LIB')
sys.path.append(LIBpath)
from libseisGT import inventory2traceid, get_FDSN_inventory, get_FDSN_Stream, removeInstrumentResponse
from libgeoGT import km2deg
from metrics import eventStatistics
from libseisGT import stream_min_starttime

# constants
#TOPDIR = '/home/seisan/MiamiLakes'
testmode =  True
TOPDIR = os.path.join(os.getenv('HOME'), 'DATA', 'MiamiLakes' )
PRE_TRIGGER_SECS = 15
POST_TRIGGER_SECS = 15
LATITUDE = 25.911667
LONGITUDE = -80.325
NETWORK = "AM"
searchRadiusKm = 20
searchRadiusDeg = km2deg(searchRadiusKm)

from obspy.clients.fdsn import Client
fdsnClient = Client("RASPISHAKE")

# Read Steve's Excel file (after conversion to CSV) into a dataFrame, df
excelfile = os.path.join('MiamiBlasts.xlsx')
df = pd.read_excel(excelfile, sheet_name='All')
df.dropna(subset=['UTCdatetime'], inplace=True)
if testmode:
    df = df.tail(4)
eventTime = list()
for index, row in df.iterrows():
    thisEventTime = obspy.core.utcdatetime.UTCDateTime(row['UTCdatetime'])
    eventTime.append(thisEventTime)
df['time'] = pd.Series(eventTime, index=df.index)
#df = df[["event_num", "time", "PPV_max", "Pa", "type", "AW-GW", "distance", "reduced_pressure", "comment"]]
df = df[["event_num", "time", "PPV_max", "Pa", "type" ]]
print(df)


# Save this improved event CSV file
csvmasterfile = os.path.join(TOPDIR, 'MiamiBlastsMasterListPandas.csv')
df.to_csv(csvmasterfile, index=False)
# # Read the CSV file that has the time column added
# df = pd.read_csv(csvmasterfile)  

def stream2wavfilename(st):
    min_stime, max_stime, min_etime, max_etime = stream_min_starttime(st)
    isofmt = min_stime.isoformat()
    wavname = isofmt[0:10] + "-" + isofmt[11:13] + isofmt[14:16] + "-" + isofmt[17:19] + 'S.MILAK__%03d' % (len(st)+1)    
    return wavname

''' ------------ Main loop over events starts here --------------'''
for index, row in df.iterrows():
    yyyymmddhhmm = row['time'].strftime('%Y%m%d%H%M')
    stationXmlFile = os.path.join(TOPDIR, 'events', 'inventories', 'Inventory.%s.xml' % yyyymmddhhmm)
    inv = get_FDSN_inventory(fdsnClient, row['time'], stationXMLfile, NETWORK, LATITUDE, LONGITUDE, searchRadiusDeg, PRE_TRIGGER_SECS, POST_TRIGGER_SECS )
            
    #inv.plot(projection='local', resolution='l')
    eventTime = list()

    # get list of unique trace ids
    trace_ids = inventory2traceid(inv)    

    # file paths to save
    rawfile = os.path.join(TOPDIR, 'events', 'mseed', 'raw', 'Event.%s.RAW.mseed' % yyyymmddhhmm)

    # load raw data and save to rawfile
    stR = get_FDSN_Stream(fdsnClient, trace_ids, rawfile, startt, endt )
    wavnameR = stream2wavfilename(stR) 
    stR.write(os.path.join(TOPDIR, 'events', 'mseed', 'detected', wavnameR), format='MSEED' ) 
    
    # detect events inside this Stream
    stD = libseisGT.detectEvent(st)
    
    # Load corrected data, or correct from raw data
    correctedfile = os.path.join(TOPDIR, 'events', 'mseed', 'corrected', wavfileR)
    if os.path.exists(correctedfile):
        # Load corrected data from file
        stC = obspy.core.read(correctedfile)
    else:
        # Correct raw data & save to file
        stC = removeInstrumentResponse(stR, preFilter = (1, 1.5, 30.0, 45.0), outputType = "VEL")
        
    # detect events inside this Stream
    stD = stC.copy()
    trig, ontimes, offtimes = libseisGT.detectEvent(stD)
    min_stime, max_stime, min_etime, max_etime = stream_min_starttime(stD)
    detectionStart = min_stime
    detectionEnd = max_etime
    if len(ontimes)>0:
        detectionStart = ontimes[0] - PRE_TRIGGER_SECS
        if detectionStart < min_stime:
            detectionStart = min_stime
    if len(offtimes>0):
        detectionEnd = offtimes[-1] + POST_TRIGGER_SECS
        if detectionEnd > max_etime:
            detectionEnd = max_etime 
    stD.trim(starttime = detectionStart, endtime = detectionEnd)
    wavnameD = stream2wavfilename(stD) 
    stD.write(os.path.join(TOPDIR, 'events', 'mseed', 'detected', wavnameD), format='MSEED' )    
    
    # Plot corrected event
    pngfile = os.path.join(TOPDIR, 'events', 'plots', 'detected', wavnameR + '.png')
    stC.plot(outfile = pngfile, equal_scale = False)
    
    # Make some measurements & save to CSV file
    statsfile = os.path.join(TOPDIR, 'events', 'statistics', '%s.csv' % wavnameD)
    if not os.path.exists(statsfile):
        eventDF = eventStatistics(stD)
        eventDF['wavfileR'] = wavnameR
        eventDF['wavfileD'] = wavnameD
        eventDF.to_csv(statsfile, index=False)
        
# Update AllEventStats.csv        
STATSDIR = "eventStats"
allDF = pd.DataFrame()
allEventFiles = sorted(glob.glob(os.path.join(STATSDIR, '*MILAK*.csv') ))
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
