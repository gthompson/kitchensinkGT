#!/usr/bin/env python
import os
import sys
from glob import glob
import numpy as np
import pandas as pd
import datetime as dt
from pprint import pprint
import matplotlib.pyplot as plt
from obspy import read, read_inventory, Stream
#from obspy.io.xseed.core import _read_resp
#from obspy.imaging.cm import obspy_sequential

LIBpath = os.path.join( os.getenv('HOME'),'src','kitchensinkGT', 'LIB')
sys.path.append(LIBpath)
from libMVO import fix_trace_id, inventory_fix_id_mvo, load_mvo_inventory \
    read_monty_wavfile_and_correct_traceIDs, enhance_stream, save_enhanced_stream, \
    metrics2df, read_enhanced_stream, plot_station_amplitude_map, parse_STATION0HYP, add_station_locations
from metrics import process_trace, choose_best_traces, select_by_index_list, ampengfft \
    Mlrichter, Eseismic_Boatwright, Eseismic2magnitude, compute_stationEnergy
from libseisGT import Stream_min_starttime, detect_network_event #, plot_seismograms 
from seisan_classes import spath2datetime, Sfile #, printEvents
from obspy.geodetics import locations2degrees, degrees2kilometers


sys.path.append(os.path.join( os.getenv('HOME'),'src', 'icewebPy') )
import IceWeb


def get_sfile_list(SEISAN_DATA, DB, startdate, enddate): 
    """
    make a list of Sfiles between 2 dates
    """

    event_list=[]
    reapath = os.path.join(SEISAN_DATA, 'REA', DB)
    years=list(range(startdate.year,enddate.year+1))
    for year in years:
        if year==enddate.year and year==startdate.year:
            months=list(range(startdate.month,enddate.month+1))
        elif year==startdate.year:
            months=list(range(startdate.month,13))
        elif year==enddate.year:
            months=list(range(1,enddate.month+1))
        else:
            months=list(range(1,13))
        for month in months:
            #print month
            yearmonthdir=os.path.join(reapath, "%04d" % year, "%02d" % month)
            flist=sorted(glob(os.path.join(yearmonthdir,"*L.S*")))
            for f in flist:
                #fdt = sfilename2datetime(f)
                fdt = spath2datetime(f)
                #print(f, fdt)
                if fdt>=startdate and fdt<enddate:
                    event_list.append(f)
    return event_list 



def tracedf2eventdf(tracedf, st, paths):
    # Summarize event
    print('Create a summary row for whole event')
    numOfRows = tracedf.shape[0]
    if bool_correct_data:
        df = tracedf[tracedf["units"] == 'm/s']
        if len(df.index)==0:
            df = tracedf
    else:
        df = tracedf
    df.sort_values(by=['quality'], inplace=True)
    df = df.head(10) # get median of 10 best rows    
    wavrow = df.median(axis = 0, skipna = True).to_dict()        
    wavrow['path']=paths['wavfile']
    wavrow['num_traces']=numOfRows
    filetime=df.iloc[0]['starttime']
    wavrow['filetime']=filetime
    try:
        wavrow['year']=filetime[0:4]
        wavrow['month']=filetime[5:7]
        wavrow['day']=filetime[8:10]
        wavrow['hour']=filetime[11:13]
        wavrow['minute']=filetime[14:16]
        wavrow['second']=filetime[17:23]
    except:
        wavrow['year']=filetime.year
        wavrow['month']=filetime.month
        wavrow['day']=filetime.day
        wavrow['hour']=filetime.hour
        wavrow['minute']=filetime.minute
        wavrow['second']=filetime.second
    
    if bool_detect_event:   
        trig, ontimes, offtimes = detect_network_event(st, sta=0.4, lta=5.0, threshon=4.0, threshoff=0.24, pad=5.0)
        print('%s: %d events detected' % (paths['picklefile'], len(ontimes)))
        durations = [t['duration'] for t in trig]
        if len(durations)>0:
            bestevent = np.argmax(durations)
            thistrig=trig[int(np.argmax(durations))]
            wavrow['ontime'] = thistrig['time']
            wavrow['offtime']=thistrig['time']+thistrig['duration']  
            wavrow['trigger_duration']=thistrig['duration']
            for item in ['coincidence_sum', 'cft_peak_wmean', 'cft_std_wmean']:
                wavrow[item]=thistrig[item]
            wavrow['detection_quality']=thistrig['coincidence_sum']*thistrig['cft_peak_wmean']*thistrig['cft_std_wmean']
            
        if bool_make_png_files:
            chosen = choose_best_traces(st, MAX_TRACES=1, include_seismic=True, 
                                include_infrasound=False, include_uncorrected=False)
            tr = st[chosen[0]] 
            plt.figure()
            plt.plot(tr.times(), tr.data)
            plt.ylabel(tr.id)   
            t0 = tr.stats.starttime
            bottom, top = plt.ylim()
            plt.vlines([wavrow['ontime']-t0, wavrow['offtime']-t0], bottom, top )
            plt.savefig(paths['detectionfile'])            

    return wavrow

def wavfile2paths(wavfile):
    paths={}
    paths['WAVDIR'] = os.path.dirname(wavfile)
    paths['wavfile'] = wavfile
    paths['wavbase'] = os.path.basename(wavfile)
    parts = paths['WAVDIR'].split('WAV')
    paths['CALDIR'] = os.path.join(parts[0],'CAL')    
    paths['miniseeddir_r'] = paths['WAVDIR'].replace('WAV', 'miniseed_r') # raw
    paths['miniseeddir_u'] = paths['WAVDIR'].replace('WAV', 'miniseed_u') # uncorrected
    paths['miniseeddir_c'] = paths['WAVDIR'].replace('WAV', 'miniseed_c') # corrected
    paths['miniseedfile_r'] = os.path.join(paths['miniseeddir_r'], paths['wavbase'] + '.mseed')
    paths['miniseedfile_u'] = os.path.join(paths['miniseeddir_u'], paths['wavbase'] + '.mseed')
    paths['miniseedfile_c'] = os.path.join(paths['miniseeddir_c'], paths['wavbase'] + '.mseed') 
    paths['traceCSVfile_u'] = paths['miniseedfile_u'].replace('.mseed', '.csv')
    paths['traceCSVfile_c'] = paths['miniseedfile_c'].replace('.mseed', '.csv')
    paths['detectionfile'] = os.path.join(paths['HTMLDIR'], paths['wavbase'] + '_detection.png')
    return paths
    

def enhanceWAV(wavfile):
    paths = wavfile2paths(wavfile)
    for item in ['miniseeddir_r', 'miniseeddir_u', 'miniseeddir_r']:
        if not os.path.exists(paths[item]):
            os.makedirs(paths[item])  

    wavrow={}
    
    ### Part I added
    wavbase = os.path.basename(wavfile)

    raw_st = read_monty_wavfile_and_correct_traceIDs(wavfile, bool_ASN=False)
    add_station_locations(raw_st, station_locationsDF)

    uncorrected_st = enhance_stream(raw_st)
    uncorrected_df = metrics2df(uncorrected_st)
    save_enhanced_stream(uncorrected_st, uncorrected_df, paths['miniseedfile_u'], save_pickle=False) # SCAFFOLD

    corrected_st = enhance_stream(raw_st, paths['CALDIR'])
    corrected_df = metrics2df(corrected_st)
    save_enhanced_stream(corrected_st, corrected_df, paths['miniseedfile_c'], save_pickle=False)
    #plot_station_amplitude_map(corrected_st, station0hypfile=station0hypfile)
    #plot_seismograms(corrected_st, outfile='corrected_seismograms.png')
    
    dome_lat = 16.7166638
    dome_lon = -62.1833326 

    mag_df = corrected_df
    mag_df['R'] = degrees2kilometers(locations2degrees(mag_df['lat'], mag_df['lon'], dome_lat, dome_lon))*1000.0
    mag_df['magA'] = Mlrichter(mag_df['peakamp'], R=mag_df['R']) 
    mag_df['Eseismic'] = Eseismic_Boatwright(mag_df['energy'], R=mag_df['R'])
    mag_df['magE'] = Eseismic2magnitude(mag_df['Eseismic'])
    for i,row in mag_df.iterrows():
        if row['units']!='m/s':
            mag_df.loc[i,'magA'] = None
            mag_df.loc[i,'Eseismic'] = None
            mag_df.loc[i,'magE'] = None
    print(mag_df)
    
    # Need to then add the station magnitudes back into the corrected_df and save it.
    mag_df.to_csv(paths['traceCSVfile_c'], index=False)

    # These stats can be used for network magnitudes
    mag_df[['magA','magE']].describe()
    
    ### End of part added

    ### Need to refactor parts below still

    st.select(component='[ENZ]') # subset to seismic components only
    if len(st)==0:
        return wavrow
   

    wavrow = tracedf2eventdf(tracedf, st, paths)        
    return wavrow    
    


def processSeisanYearMonth(SEISAN_DATA, SEISAN_DB, YYYY, MM, filesdone, MAXFILES=999999, startdate=None):
    failedWAVfiles=[]
    LoD = []
    
    # We aim to add a couple of columns from the S-file, and save to this
    reawavCSVfile=os.path.join(SEISAN_DATA, 'reawav_%s%s%s.csv' % (SEISAN_DB, YYYY, MM) )
    if os.path.exists(reawavCSVfile):
        return failedWAVfiles, filesdone
    
    # Get s-file list
    if not startdate:
        startdate = dt.datetime(int(YYYY), int(MM), 1)
    if int(MM)<12:
        enddate = dt.datetime(int(YYYY), int(MM)+1, 1)
    else:
        enddate = dt.datetime(int(YYYY)+1, 1, 1)
    print(startdate, enddate)
    return
    slist = sorted(get_sfile_list(SEISAN_DATA, SEISAN_DB, startdate, enddate))
    
    for i,sfile in enumerate(slist):
        print('Processing %d of %d: %s' % (i, len(slist), sfile) )
        if i==MAXFILES:
            break
            
        s = Sfile(sfile, use_mvo_parser=True)
        #s.cat()
        #s.printEvents()
        d = s.to_dict()
        #pprint(d)
        
        for item in ['wavfile1', 'wavfile2']:
            if d[item]:
                if os.path.exists(d[item]):
                    wavbase = os.path.basename(d[item])
                    if 'MVO' in wavbase:
                        print('Processing ',d[item])
                        eventrow=[]
                        eventrow = enhanceWAV(d[item])
                        if eventrow:
                            eventrow['sfile']=os.path.basename(s.path)
                            eventrow['mainclass']=s.mainclass
                            eventrow['subclass']=s.subclass
                            LoD.append(eventrow)
                            filesdone += 1
                            if filesdone >= MAXFILES:
                                break
    if LoD:
        df = pd.DataFrame(LoD)
        print('Writing ',reawavCSVfile)
        df.drop(df.filter(regex="Unname"),axis=1, inplace=True)
        df = df.set_index('filetime')       
        df.to_csv(reawavCSVfile, index=True) 
    return failedWAVfiles, filesdone

#######################################################################

SEISAN_DATA = os.path.join( os.getenv('HOME'),'DATA','MVO')
os.chdir(SEISAN_DATA)
SEISAN_DB = 'MVOE_'
station0hypfile = os.path.join(SEISAN_DATA, 'DAT', 'STATION0_MVO.HYP')
station_locationsDF = parse_STATION0HYP(station0hypfile) 
bool_shortperiod=False
bool_correct_data=True
bool_make_png_files=False
bool_detect_event=True
bool_overwrite=False
#startdate = dt.datetime(2007,1,29,0,0,0)
startdate = dt.datetime(1996,10,20,0,0)

filesdone = 0
yeardirs = sorted(glob(os.path.join('REA',SEISAN_DB,'[12]???')))
for yeardir in yeardirs:
    YYYY = os.path.basename(yeardir)
    monthsdirs = sorted(glob(os.path.join(yeardir,'[01]?')))
    for monthdir in monthsdirs:
        if filesdone>=MAXFILES:
            break
        MM = os.path.basename(monthdir)
        if startdate:
            SYYYY = startdate.year
            SMM = startdate.month                
            if SYYYY>int(YYYY):
                continue
            elif SMM > int(MM):
                continue
        print('\n**** Processing %s ****' % monthdir)
        failedWAVfiles, filesdone = processSeisanYearMonth('.', SEISAN_DB, YYYY, MM, filesdone, MAXFILES=999999, startdate=startdate)
        if len(failedWAVfiles)>0:
            fptr=open('failedWAVfiles.txt','a')
            for element in failedWAVfiles:
                fptr.write(element)
                fptr.write('\n')
            fptr.close() 
