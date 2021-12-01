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
from metrics import process_trace, select_by_index_list, ampengfft, \
    Mlrichter, Eseismic_Boatwright, Eseismic2magnitude, compute_stationEnergy
from libseisGT import Stream_min_starttime #, plot_seismograms 
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
    

def enhanceWAV(wavfile, bool_overwrite=False, station_locationsDF=None, MASTER_INV=None):
    bool_do_everything = bool_overwrite
    paths = wavfile2paths(wavfile)
    for item in ['miniseeddir_r', 'miniseeddir_u', 'miniseeddir_r']:
        if not os.path.exists(paths[item]):
            os.makedirs(paths[item]) 

    wavbase = os.path.basename(wavfile)

    
    # we only need to create the uncorrected miniseed file if it does not exist - which means we need to load the wavfile first
    if (not os.path.exists(paths['miniseedfile_u'])) or bool_overwrite:
        bool_do_everything = True
        print('Loading %s' % wavfile)
    
        raw_st = read_monty_wavfile_and_correct_traceIDs(wavfile, bool_ASN=False)
        add_station_locations(raw_st, station_locationsDF)

        uncorrected_st = enhance_stream(raw_st)

        uncorrected_df = metrics2df(uncorrected_st)
        save_enhanced_stream(uncorrected_st, uncorrected_df, paths['miniseedfile_u'], save_pickle=False) # SCAFFOLD

    

    # we only need to create the corrected miniseed file if it does not exist - which means we need to load the wavfile first
    # unless we did that above (does raw_st exist?)
    if (not os.path.exists(paths['miniseedfile_c'])) or bool_overwrite: 
        print('Creating %s' % paths['miniseedfile_c'])
        try:
            raw_st
        except:
            print('Loading %s' % (wavfile))  
            raw_st = read_monty_wavfile_and_correct_traceIDs(wavfile, bool_ASN=False)
            add_station_locations(raw_st, station_locationsDF)
        
        corrected_st = enhance_stream(raw_st, paths['CALDIR'], master_inv=MASTER_INV)
    
    if bool_do_everything:
        # which traces are in raw_st but not corrected_st, and why?
        corrected_ids = [tr.id for tr in corrected_st]
        uncorrected_ids = [tr.id for tr in uncorrected_st]
        for tr in raw_st:
            if tr.id in uncorrected_ids:
                if not tr.id in corrected_ids:
                    print('%s: could not correct %s' % (wavfile, tr.id))
            else:
                print('%s: poor trace %s' % (wavfile, tr.id))
                

    # we only need to create the corrected CSV file if it does not exist 
    # we may have to load the corrected miniseed file, unless we did that above
    if (not os.path.exists(paths['traceCSVfile_c'])) or bool_overwrite:
        print('Creating %s' % paths['traceCSVfile_c'])        
        try:
            corrected_st
        except:
            if os.path.exists(paths['miniseedfile_c']):
                corrected_st = read(paths['miniseedfile_c'])

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
        #mag_df[['magA','magE']].describe()
        
    return paths['miniseedfile_c']
    
def process1event(sfile, bool_overwrite, station_locationsDF=None,  MASTER_INV=None, bool_index_only=False):
    
    try:
        s = Sfile(sfile, fast_mode=True)
        d = s.to_dict()
        sfileindex_dict = {'sfile':os.path.basename(s.path), 'DSN_wavfile':None, 'DSN_exists':False, 'ASN_wavfile':None, 'ASN_exists':False, 'corrected_DSN_mseed':None, 'corrected_ASN_mseed':None, 'mainclass':s.mainclass, 'subclass':s.subclass}
    except:
        os.system('echo %s >> bad_sfiles.log' % sfile)
        return None
    #s.cat()
    #s.printEvents()        
    #pprint(d) 
         
    for item in ['wavfile1', 'wavfile2']:
        if d[item]:
            wavbase = os.path.basename(d[item])
            if 'MVO' in wavbase:
                sfileindex_dict['DSN_wavfile'] = wavbase
                if os.path.exists(d[item]):
                    sfileindex_dict['DSN_exists'] = True
                    if bool_index_only:
                        paths = wavfile2paths(d[item])  
                        DSN_mseedfile = paths['miniseedfile_c']
                        if os.path.exists(DSN_mseedfile):
                            sfileindex_dict['corrected_DSN_mseed'] = DSN_mseedfile
                    else:
                        try:
                            DSN_mseedfile = enhanceWAV(d[item], bool_overwrite=bool_overwrite, station_locationsDF=station_locationsDF,  MASTER_INV=MASTER_INV)
                            sfileindex_dict['corrected_DSN_mseed'] = DSN_mseedfile
                        except:
                            os.system('echo %s, %s >> seisan2pandas_failed.log' % (sfile, wavbase))

            elif 'SPN' in wavbase:
                sfileindex_dict['ASN_wavfile'] = wavbase
                if os.path.exists(d[item]):
                    sfileindex_dict['ASN_exists'] = True 
    return sfileindex_dict
        

def set_globals():
    SEISAN_DATA = os.path.join( os.getenv('HOME'),'DATA','MVO')
    #master_station_xml = './MontserratDigitalSeismicNetwork.xml'
    #os.system("cp %s %s/CAL/" % (master_station_xml, SEISAN_DATA) )
    os.chdir(SEISAN_DATA)
    master_station_xml = 'CAL/MontserratDigitalSeismicNetwork.xml'
    SEISAN_DB = 'MVOE_'
    station0hypfile = os.path.join(SEISAN_DATA, 'DAT', 'STATION0_MVO.HYP')
    if os.path.exists(station0hypfile):
        station_locationsDF = parse_STATION0HYP(station0hypfile) 
    else:
        station_locationsDF = None

    MASTER_INV = read_inventory(master_station_xml)
    bool_overwrite=False
    return SEISAN_DATA, SEISAN_DB, station_locationsDF, MASTER_INV, bool_overwrite

