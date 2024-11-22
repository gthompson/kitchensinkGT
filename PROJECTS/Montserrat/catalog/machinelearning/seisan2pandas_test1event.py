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
from libMVO import fix_trace_id, inventory_fix_id_mvo, load_mvo_inventory, \
    read_monty_wavfile_and_correct_traceIDs, enhance_stream, save_enhanced_stream, \
    metrics2df, read_enhanced_stream, plot_station_amplitude_map, parse_STATION0HYP, \
    add_station_locations, load_mvo_master_inventory
from metrics import process_trace, choose_best_traces, select_by_index_list, ampengfft, \
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



def eventdf2catalogdf(eventdf, st, paths):
    # Summarize event
    print('Create a summary row for whole event')
    numOfRows = eventdf.shape[0]
    if bool_correct_data:
        df = eventdf[eventdf["units"] == 'm/s']
        if len(df.index)==0:
            df = eventdf
    else:
        df = eventdf
    df.sort_values(by=['quality'], inplace=True)
    df = df.head(10) # get 10 best rows    
    catalogrow = df.median(axis = 0, skipna = True).to_dict()        
    catalogrow['path']=paths['wavfile']
    catalogrow['num_traces']=numOfRows
    filetime=df.iloc[0]['starttime']
    catalogrow['filetime']=filetime
    try: # is filetime a string or a datetime object?
        catalogrow['year']=filetime[0:4]
        catalogrow['month']=filetime[5:7]
        catalogrow['day']=filetime[8:10]
        catalogrow['hour']=filetime[11:13]
        catalogrow['minute']=filetime[14:16]
        catalogrow['second']=filetime[17:23]
    except:
        catalogrow['year']=filetime.year
        catalogrow['month']=filetime.month
        catalogrow['day']=filetime.day
        catalogrow['hour']=filetime.hour
        catalogrow['minute']=filetime.minute
        catalogrow['second']=filetime.second
    
    if bool_detect_event:   
        trig, ontimes, offtimes = detect_network_event(st, sta=0.4, lta=5.0, threshon=4.0, threshoff=0.24, pad=5.0)
        print('%d events detected' % (len(ontimes)))
        durations = [t['duration'] for t in trig]
        if len(durations)>0:
            bestevent = np.argmax(durations) # choose longest detection
            thistrig=trig[int(np.argmax(durations))]
            catalogrow['ontime'] = thistrig['time']
            catalogrow['offtime']=thistrig['time']+thistrig['duration']  
            catalogrow['trigger_duration']=thistrig['duration']
            for item in ['coincidence_sum', 'cft_peak_wmean', 'cft_std_wmean']:
                catalogrow[item]=thistrig[item]
            catalogrow['detection_quality']=thistrig['coincidence_sum']*thistrig['cft_peak_wmean']*thistrig['cft_std_wmean']
            
            chosen = choose_best_traces(st, MAX_TRACES=1, include_seismic=True, 
                                include_infrasound=False, include_uncorrected=False)
            tr = st[chosen[0]] 
            plt.figure()
            plt.plot(tr.times(), tr.data)
            plt.ylabel(tr.id)   
            t0 = tr.stats.starttime
            bottom, top = plt.ylim()
            plt.vlines([catalogrow['ontime']-t0, catalogrow['offtime']-t0], bottom, top )
            plt.savefig(paths['detectionfile'])            

    return catalogrow

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
    paths['detectionfile'] = os.path.join(paths['miniseeddir_c'], paths['wavbase'] + '_detection.png')
    return paths
    

def enhanceWAV(wavfile):
    paths = wavfile2paths(wavfile)
    for item in ['miniseeddir_r', 'miniseeddir_u', 'miniseeddir_r']:
        if not os.path.exists(paths[item]):
            os.makedirs(paths[item])  
    
    ### Part I added
    wavbase = os.path.basename(wavfile)

    raw_st = read_monty_wavfile_and_correct_traceIDs(wavfile, bool_ASN=False)
    add_station_locations(raw_st, station_locationsDF)

    uncorrected_st = enhance_stream(raw_st)
    uncorrected_df = metrics2df(uncorrected_st)
    save_enhanced_stream(uncorrected_st, uncorrected_df, paths['miniseedfile_u'], save_pickle=False) # SCAFFOLD

    corrected_st = enhance_stream(raw_st, paths['CALDIR'], master_inv=MASTER_INV)
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

    corrected_st.select(component='[ENZ]') # subset to seismic components only
    if len(corrected_st)==0:
        eventrow = {}
        return eventrow
    else:
        eventrow = eventdf2catalogdf(mag_df, corrected_st, paths)        
        return eventrow    
    


def test1event(sfile):
          
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


#######################################################################

SEISAN_DATA = os.path.join( os.getenv('HOME'),'DATA','MVO')
master_station_xml = './MontserratDigitalSeismicNetwork.xml'
os.system("cp %s %s/CAL/" % (master_station_xml, SEISAN_DATA) )
master_station_xml = 'CAL/MontserratDigitalSeismicNetwork.xml'
os.chdir(SEISAN_DATA)
SEISAN_DB = 'MVOE_'
station0hypfile = os.path.join(SEISAN_DATA, 'DAT', 'STATION0_MVO.HYP')
if os.path.exists(station0hypfile):
    station_locationsDF = parse_STATION0HYP(station0hypfile) 
else:
    print(station0hypfile,' not found')
    exit()
print(station_locationsDF)
MASTER_INV = read_inventory(master_station_xml)
bool_shortperiod=False
bool_correct_data=True
bool_make_png_files=False
bool_detect_event=True
bool_overwrite=False
sfilebase = '04-0556-28L.S199703'
test1event(os.path.join(SEISAN_DATA,'REA',SEISAN_DB,sfilebase[-6:-2],sfilebase[-2:],sfilebase))

