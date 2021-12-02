#!/usr/bin/env python
import os
import sys
from glob import glob
import numpy as np
import pandas as pd
import datetime as dt
from obspy import Stream  
LIBpath = os.path.join( os.getenv('HOME'),'src','kitchensinkGT', 'LIB')
sys.path.append(LIBpath)
from libSeisan2Pandas import set_globals

from libseisGT import detect_network_event
from metrics import choose_best_traces


def create_catalog(SEISAN_DATA, SEISAN_DB, YYYY, MM, bool_detect_event):
    # the details file created here is what would be used for fingerprints, and forms the actual catalog.
    # this would however be much faster to generate if we were not trying to detect events, because
    # then there would be no need to read the miniseed file and run the detection
    miniseeddir = os.path.join(SEISAN_DATA, 'miniseed_c', SEISAN_DB)
    sfileindexfile=os.path.join(miniseeddir, 'sfileindex_%s%s%s.csv' % (SEISAN_DB, YYYY, MM) )
    if not os.path.exists(sfileindexfile):
        print(sfileindexfile, ' does not exist')
        return
    print(sfileindexfile, ' found')
    sfileindex_df = pd.read_csv(sfileindexfile)
    print('sfile index for %s/%s has %d rows' % (YYYY,MM,len(sfileindex_df)))
    if bool_detect_event:
        catalogfile=os.path.join(miniseeddir, 'detected_catalog_%s%s%s.csv' % (SEISAN_DB, YYYY, MM) )
    else:
        catalogfile=os.path.join(miniseeddir, 'catalog_%s%s%s.csv' % (SEISAN_DB, YYYY, MM) )
        for i,row in sfileindex_df.iterrows():
            #print(row)
            mseedfile = row['corrected_DSN_mseed']
            if not isinstance(mseedfile, str):
                continue
            if not os.path.exists(mseedfile): # corrected_DSN_mseed is full path, so if sfileindex came from a different computer, start of path may be wrong
                fileparts = mseedfile.split(SEISAN_DB)
                #print(fileparts)
                parentdir = os.path.basename(fileparts[0][:-1])
                #print(SEISAN_DATA, parentdir, SEISAN_DB)
                mseedfile = os.path.join(SEISAN_DATA, parentdir, SEISAN_DB) +  fileparts[1]
                #print(mseedfile)
            
            # this part seems to be same functionality as read_enhanced_wav
            if isinstance(mseedfile, str):
                
                eventcsvfile = mseedfile.replace('.mseed','.csv')
                event_df = pd.read_csv(eventcsvfile)
                if len(event_df.index)==0:
                    continue
                
                corrected_st = Stream()
                if bool_detect_event:
                    corrected_st = read(mseedfile)
                    corrected_st.select(component='[ENZ]') # subset to seismic components only    

                catalogrow = eventdf2catalogrow(event_df, mseedfile, bool_detect_event)
                        
                # add items from the event summary to the index summary
                for key in catalogrow.keys():
                    sfileindex_df.loc[i,key] = catalogrow[key]
                 
        sfileindex_df.to_csv(catalogfile)          

      
def eventdf2catalogrow(df, mseedfile, bool_detect_event):
    # Summarize event
    #print('Create a summary row for whole event')
    df.sort_values(by=['quality'], inplace=True)
    df = df.head(10) # get 10 best rows    
    catalogrow = df.median(axis = 0, skipna = True).to_dict()        
    catalogrow['path']=os.path.basename(mseedfile).replace('.mseed','')
    catalogrow['num_traces']=len(df.index)
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
    
    # Optionally add detection metrics and make a plot - but time consuming
    if bool_detect_event and isinstance(st, Stream):   
        corrected_st = read(mseedfile)
        corrected_st.select(component='[ENZ]') # subset to seismic components only
        if len(corrected_st) > 0:
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
                detectionpngfile = waveformfilepath.replace('.mseed','_detection.png')
                plt.savefig(detectionpngfile)            

    return catalogrow        

#######################################################################
print('Usage %s [STARTYEAR [STARTMONTH [MAXFILES] ] ]' % sys.argv[0])
STARTYEAR = "1996"
STARTMONTH = "01"
MAXFILES = 999999
argc = len(sys.argv)
print(argc)
if argc>1:
    STARTYEAR = sys.argv[1]
if argc>2:
    STARTMONTH = sys.argv[2]
if argc>3:
    MAXFILES = sys.argv[3]
    
# import controlling variables - actually not globals
SEISAN_DATA, SEISAN_DB, station_locationsDF, MASTER_INV, bool_overwrite = set_globals()
if not isinstance(station_locationsDF, pd.DataFrame):
    print('Station coordinates not loaded. Exiting. Although we could change code and get these from StationXML file')
    exit()
    
# the default is to loop over years/months
filesdone = 0
bool_detect_event = False
yeardirs = sorted(glob(os.path.join(SEISAN_DATA, 'REA',SEISAN_DB,'[12]???')))
for yeardir in yeardirs:
    print(yeardir)
    YYYY = os.path.basename(yeardir)
    if STARTYEAR:
        if YYYY<STARTYEAR:
            continue
    monthsdirs = sorted(glob(os.path.join(yeardir,'[01]?')))
    for monthdir in monthsdirs:
        if filesdone>=MAXFILES:
            break
        MM = os.path.basename(monthdir)
        if STARTYEAR and STARTMONTH>"01":
            if YYYY==STARTYEAR and MM<STARTMONTH:
                continue       
        print('\n**** Processing %s ****' % monthdir)
        create_catalog(SEISAN_DATA, SEISAN_DB, YYYY, MM, bool_detect_event)

