#!/usr/bin/env python
# coding: utf-8

# # Segment SDS Archive based on launch event times and generate event web browser
# 
# Launch times come from 'PilotStudy_KSC_Rocket_Launches.xlsx'
# SDS archives are at SDS_TOP and contain data from wells 6I and 6S and from seismo-acoustic stations
# Segmented event waveform files are saved as MiniSEED to EVENT_WAVEFORMS
# 
# 
#from IPython import get_ipython
#get_ipython().run_line_magic('run', 'header.ipynb')
import os
import header
paths = header.setup_environment()
paths['SDS_TOP'] = os.path.join(paths['outdir'], 'SDS')
HTML_DIR = '/var/www/html/thompsong/KSC_EROSION/EVENTS'
PNG_DIR = os.path.join(HTML_DIR, 'images')
EVENT_WAVEFORMS = os.path.join(paths['outdir'], 'EVENTS') # must exist, and Excel file must be here
csv_launches = os.path.join(paths['outdir'], 'PilotStudy_KSC_Rocket_Launches.csv')
csv_launches_detected = os.path.join(paths['outdir'], 'PilotStudy_KSC_Rocket_Launches_detected.csv')
#import glob
import numpy as np
import pandas as pd
from obspy.core import read, Stream, UTCDateTime
#import FDSNtools
#import wrappers
import SDS
import libWellData as LLE

import libDatascopeGT
from obspy.clients.filesystem.sds import Client
def try_different_waveform_loading_methods(thisSDSobj, startt, endt, methods=[1,2,3,4]):

    # Datascope
    dbpath = os.path.join(paths['outdir'], 'db', f'db{startt.strftime("%Y%m%d")}' )
    #st1 = libDatascopeGT.db2Stream(dbpath, startt, endt)

    # Direct Miniseed load with ObsPy
    st2 = libDatascopeGT.mseed2Stream(paths['SDS_TOP'], startt, endt)

    # ObsPy SDS archive reader
    sdsclient = Client(paths['SDS_TOP'])
    st3 = sdsclient.get_waveforms("*", "*", "*", "[HDCES]*", startt, endt)

    # My SDS class that wraps ObsPy SDS reader
    thisSDSobj.read(startt, endt, speed=1)
    st4 = thisSDSobj.stream
      
    #compare_trace_ids(st1, st2, testA=1, testB=2)
    #compare_trace_ids(st3, st2, testA=3, testB=2)
    #compare_trace_ids(st4, st2, testA=4, testB=2) 
    #combine_streams(st2, st1)
    combine_streams(st2, st3)    
    combine_streams(st2, st4)
    return st2
    
def combine_streams(stB, stA):
    appended = False
    for trA in stA:
        found = False
        for trB in stB:
            if trA.stats.station == trB.stats.station and trA.stats.location == trB.stats.location and trA.stats.channel == trB.stats.channel:
                if trA.stats.network == '':
                    trA.stats.network = trB.stats.network
                if trA.stats.starttime >= trB.stats.starttime and trA.stats.endtime <= trB.stats.endtime:
                    found = True
                    break
        if not found:
            stB.append(trA)
            appended = True
    if appended:
        stB.merge(method=0, fill_value=0)



def event_sds2pkl(launchtime, thisSDSobj, EVENT_WAVEFORMS, pretrig=3600, posttrig=3600, overwrite=False):    
    rawfile = os.path.join(EVENT_WAVEFORMS, '%s_raw.pkl' % launchtime.strftime('%Y%m%dT%H%M%S'))
    if os.path.exists(rawfile) and not overwrite:
        print('%s already exists' % rawfile)
    else:
        print('segmenting %s from SDS' % rawfile) 
        startt = launchtime - pretrig
        endt = launchtime + posttrig
        st = try_different_waveform_loading_methods(thisSDSobj, startt, endt, methods=[1,2,3,4])

        if len(st)>0:
            try:
                st.write(rawfile, format='pickle')
            except:
                st2 = Stream()
                for tr in st:
                    try:
                        tr.write('tmp.pkl', 'pickle')
                        st2.append(tr)
                    except:
                        print('Failed:\n',tr)
                st=st2
                if len(st)>0:
                    st.write(rawfile, format='pickle')
                else:
                    rawfile=None
        else:
            print('Got no data')
            rawfile = None
    if rawfile:
        rawpng = os.path.join(PNG_DIR, os.path.basename(rawfile.replace('.pkl','.png')))
        if not os.path.exists(rawpng):
            print('creating raw file plot ',rawpng)
            st.plot(equal_scale=False, outfile=rawpng)
    return rawfile

def clean(st):
    for tr in st:
        if tr.stats.network != 'FL':
            continue
        tr.detrend('linear')
        tr.filter('highpass', freq=0.2, corners=2) 
        
def apply_calibration_correction(st):
    # calibration correction

    for tr in st:
        if 'countsPerUnit' in tr.stats:
            continue
        else:
            tr.stats['countsPerUnit'] = 1
            if not 'units' in tr.stats:
                tr.stats['units'] = 'Counts'
            if tr.stats.station[0].isnumeric(): # well data
                if len(tr.stats.network)==0:
                    tr.stats.network = '6'
                if tr.stats.channel[2] == 'D':
                    tr.stats.countsPerUnit = 1/LLE.psi2inches(1) # counts (psi) per inch
                    tr.stats.units = 'inches'
                elif tr.stats.channel[2] == 'H':
                    tr.stats.countsPerUnit = 1/6894.76 # counts (psi) per Pa
                    tr.stats.units = 'Pa'
            elif tr.stats.channel[1]=='D':
                tr.stats.countsPerUnit = 720 # counts/Pa on 1 V FS setting
                if tr.id[:-1] == 'FL.BCHH3.10.HD':
                    if tr.stats.starttime < UTCDateTime(2022,5,26): # Chaparral M25. I had it set to 1 V FS. Should have used 40 V FS. 
                        if tr.id == 'FL.BCHH3.10.HDF':
                            tr.stats.countsPerUnit = 8e5 # counts/Pa
                        else:
                            tr.stats.countsPerUnit = 720 # counts/Pa 
                    else: # Chaparral switched to 40 V FS
                        if tr.id == 'FL.BCHH3.10.HDF':
                            tr.stats.countsPerUnit = 2e4 # counts/Pa
                        else:
                            tr.stats.countsPerUnit = 18 # counts/Pa 
                tr.stats.units = 'Pa'

            elif tr.stats.channel[1]=='H':
                tr.stats.countsPerUnit = 3e2 # counts/(um/s)
                tr.stats.units = 'um/s'
            tr.data = tr.data/tr.stats.countsPerUnit
    
def maxamp(tr):
    return np.max(np.abs(tr.data))

def remove_spikes(st):
    SEISMIC_MAX = 0.1 # m/s
    INFRASOUND_MAX = 3000 # Pa
    FEET_MAX = 21 # feet
    #SEISMIC_MIN = 1e-9
    #INFRASOUND_MIN = 0.01
    
    for tr in st:
        ma = maxamp(tr)
        if tr.stats.units == 'm/s':
            tr.data[tr.data > SEISMIC_MAX] = np.nan
            tr.data[tr.data < -1 * SEISMIC_MAX] = np.nan             
        elif tr.stats.units == 'Pa':
            tr.data[tr.data > INFRASOUND_MAX] = np.nan
            tr.data[tr.data < -1 * INFRASOUND_MAX] = np.nan   
        elif tr.stats.units == 'feet':
            tr.data[tr.data > FEET_MAX] = np.nan
            tr.data[tr.data < -1 * FEET_MAX] = np.nan               

from obspy.signal.trigger import coincidence_trigger
from pprint import pprint
import matplotlib.dates as dates
def detectEvent(st, launchtime):
    trig = coincidence_trigger("recstalta", 3.5, 1, st, 3, sta=2, lta=40)
    best_trig = {}
    best_product = 0
    for this_trig in trig:
        thistime = dates.date2num(this_trig['time'])
        this_product = this_trig['coincidence_sum']*this_trig['duration']
        if this_product > best_product:
            best_trig = this_trig
            best_product = this_product
    pprint(best_trig)
    return best_trig['time']

'''
def add_snr(st, assoctime, threshold=1.5):
    nstime = max([st[0].stats.starttime, assoctime-240])
    netime = min([st[0].stats.endtime, assoctime-60])
    sstime = assoctime
    setime = min([st[0].stats.endtime, assoctime+120])    
    for tr in st:
        tr_noise = tr.copy().trim(starttime=nstime, endtime=netime)
        tr_signal = tr.copy().trim(starttime=sstime, endtime=setime)
        tr.stats['noise'] = np.nanmedian(np.abs(tr_noise.data))
        tr.stats['signal'] = np.nanmedian(np.abs(tr_signal.data))
        tr.stats['snr'] = tr.stats['signal']/tr.stats['noise']
        '''

def group_streams_for_plotting(st):
    groups = {}
    stationsWELL = ['6S', '6I']
    for station in stationsWELL:
        stationStream = st.select(network=station)
        #stationIDS = list(set([tr.id for tr in stationStream]))
        groups[station] = stationStream
    streamSA = st.select(network='FL')
    stationsSA = list(set([tr.stats.station for tr in streamSA]))
    for station in stationsSA:
        stationStream = streamSA.select(station=station)
        #stationIDS = list(set([tr.id for tr in stationStream]))
        groups[station] = stationStream
    #print(groups)
    return groups  




# Read launch data into a DataFrame and generate a list of launch times in Python datetime.datetime format


startover = True # starts with original CSV file again
if os.path.isfile(csv_launches_detected) and startover==False:
    launchesDF = LLE.removed_unnamed_columns(pd.read_csv(csv_launches_detected, index_col=None))
else:
    launchesDF = LLE.removed_unnamed_columns(pd.read_csv(csv_launches, index_col=None))
    dt_tmp = pd.to_datetime(launchesDF['Date'] + ' ' +  launchesDF['Time'])
    launchesDF['Date'] = [pdt.to_pydatetime() for pdt in dt_tmp]
    launchesDF.drop(labels='Time', axis=1, inplace=True)
    del dt_tmp

for thisdir in [EVENT_WAVEFORMS, HTML_DIR, PNG_DIR]:
    if not os.path.isdir(thisdir):
        os.makedirs(thisdir)

if not 'rawfile' in launchesDF.columns:
    launchesDF['rawfile'] = ''
if not 'corrected_file' in launchesDF.columns: 
    launchesDF['corrected_file'] = ''
if not 'detection_time' in launchesDF.columns:
    launchesDF['detection_time'] = '' 
if not 'short_file' in launchesDF.columns:
    launchesDF['short_file'] = ''
if not 'plotted' in launchesDF.columns:
    launchesDF['plotted'] = False     
launchesDF.to_csv(csv_launches_detected) 

def check_st(st):
    print(st)
    wellst=st.select(network='6*')
    if len(wellst)>0:
        plot(equal_scale=False)
        anykey = input('<Enter> to contine')

# For each launch, segment raw SDS data to multi-trace MiniSEED file in EVENT_WAVEFORMS directory
thisSDSobj = SDS.SDSobj(paths['SDS_TOP'])
print(launchesDF)
for i, row in launchesDF.iterrows():
    launchTime = UTCDateTime(row['Date'])
    if not row['corrected_file']:
        print('Processing launch at %s' % launchTime.strftime('%Y-%m-%d %H:%M:%S')) 
        print(row)
        if row['rawfile']:
            rawfile = os.path.join(EVENT_WAVEFORMS, row['rawfile'])
        else:
            rawfile = event_sds2pkl(launchTime, thisSDSobj, EVENT_WAVEFORMS, overwrite=False)
            if rawfile:
                launchesDF.at[i, 'rawfile'] = os.path.basename(rawfile)
            else:
                raise Exception('failed to create %s' % rawfile)
        try:
            st = read(rawfile)    
            #st.merge(method=0, fill_value=0)
        except:
            st = Stream()
        print('%s: %d channels' % (rawfile,len(st)))

        # all these functions safe for well traces too
        print(st)
        st.plot
        clean(st) 
        print(st)
        apply_calibration_correction(st)
        print(st)
        remove_spikes(st)
        print(st)
        
        # write corrected event out
        correctedfile =  os.path.join(EVENT_WAVEFORMS, '%s_long.pkl' % launchTime.strftime('%Y%m%dT%H%M%S'))
        print('Writing %s' % correctedfile)
        try:
            st.write(correctedfile, format='PICKLE') # save 2-hour event waveforms
            launchesDF.at[i, 'corrected_file'] = os.path.basename(correctedfile)
        except:
            pass
del thisSDSobj
launchesDF.to_csv(csv_launches_detected) 


           
for i, row in launchesDF.iterrows():
    if row['detection_time'] or not startover:
        continue
    launchTime = UTCDateTime(row['Date'])
    if row['corrected_file']:
        correctedfile =  os.path.join(EVENT_WAVEFORMS, row['corrected_file'])
    else:
        continue
    print('Detecting launch at %s' % launchTime.strftime('%Y-%m-%d %H:%M:%S'))                        

    # subset out the seismo-acoustic traces for detection purposes
    if not os.path.isfile(correctedfile):
        print('File not found: ',correctedfile)
        continue
    st = read(correctedfile)
    SA = st.copy().select(network='FL').trim(starttime=launchTime-100, endtime=launchTime+200)
    if len(SA)==0:
        continue

    assocTime = detectEvent(SA, launchTime)
            
    if abs(assocTime-launchTime)>100:
        assocTime=launchTime
    launchesDF.at[i, 'detection_time'] = assocTime
launchesDF.to_csv(csv_launches_detected) 


# Short files
for i, row in launchesDF.iterrows():                       
    if not row['short_file']:
        launchTime = UTCDateTime(row['Date']) 
        print('- Creating short file for launch at %s' % launchTime.strftime('%Y-%m-%d %H:%M:%S')) 
        if row['corrected_file']:
            correctedfile =  os.path.join(EVENT_WAVEFORMS, row['corrected_file'])
        else:
            continue   
        if not os.path.isfile(correctedfile):
            print('File not found: ',correctedfile)
            continue
        # save 3-minute event waveforms
        st = read(correctedfile)
        if len(st)==0:
            continue

        assocTime = row['detection_time']
        if not assocTime:
            continue
        st_short = st.copy()
        #st_short.filter('highpass', freq=0.1, corners=2)
        st_short.trim(starttime=assocTime-30, endtime=assocTime+150)
        print(st_short)
        if len(st_short)>0:
            # write shorter corrected event out
            shortfile =  os.path.join(EVENT_WAVEFORMS, '%s_short.pkl' % launchTime.strftime('%Y%m%dT%H%M%S'))
            print('Writing %s' % shortfile)
            try:
                st_short.write(shortfile, format='PICKLE') # save 2-hour event waveforms
            except:
                pass
            launchesDF.at[i, 'short_file'] = os.path.basename(shortfile)
launchesDF.to_csv(csv_launches_detected)         


# Plots
for i, row in launchesDF.iterrows():
    if row['plotted']:
        continue
    launchTime = UTCDateTime(row['Date'])
    print('- Plotting launch at %s' % launchTime.strftime('%Y-%m-%d %H:%M:%S')) 
    for ext in ['long', 'short']:
        if ext=='short':
            pklfile = row['short_file']
        else:    
            pklfile = row['corrected_file']
        if not pklfile:
            continue
        pklfile = os.path.join(EVENT_WAVEFORMS, pklfile)
        if not os.path.isfile(pklfile):
            print('File not found: ',pklfile)
            continue
        st = read(pklfile)
        if len(st)==0:
            continue        
        groups = group_streams_for_plotting(st)
        for station, stream_group in groups.items():
            if len(stream_group)>0:
                pngfile = os.path.join(PNG_DIR, '%s_%s_%s.png' % (launchTime.strftime('%Y%m%dT%H%M%S'), station, ext))
                stream_group.plot(equal_scale=False, outfile=pngfile)
    launchesDF.at[i, 'plotted'] = True
launchesDF.to_csv(csv_launches_detected)          

