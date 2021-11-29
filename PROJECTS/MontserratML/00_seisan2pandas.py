#!/usr/bin/env python
import os
import sys
from glob import glob
import numpy as np
import pandas as pd
import datetime as dt
LIBpath = os.path.join( os.getenv('HOME'),'src','kitchensinkGT', 'LIB')
sys.path.append(LIBpath)
from libSeisan2Pandas import get_sfile_list, process1event, set_globals
from libseisGT import detect_network_event
from metrics import choose_best_traces

    
def processSeisanYearMonth(SEISAN_DATA, SEISAN_DB, YYYY, MM, filesdone, MAXFILES=999999, startdate=None, \
                           bool_overwrite=False, station_locationsDF=None, MASTER_INV=None):
    #summaryLoD = []
    summaryfile=os.path.join(SEISAN_DATA, 'miniseed_c', SEISAN_DB, 'summary_%s%s%s.csv' % (SEISAN_DB, YYYY, MM) )
    if os.path.exists(summaryfile) and bool_overwrite:
        summarydf = pd.read_csv(summaryfile)
    else:
        summarydf = pd.DataFrame(columns=['sfile', 'DSN_wavfile', 'DSN_exists', 'ASN_wavfile', 'ASN_exists', 'corrected_DSN_mseed', 'corrected_ASN_mseed', 'mainclass', 'subclass']) 

    # Get s-file list
    startdate_thismonth = dt.datetime(int(YYYY), int(MM), 1)
    if not startdate:
        startdate = startdate_thismonth
    else:
        if startdate < startdate_thismonth:
            startdate = startdate_thismonth
    if int(MM)<12:
        enddate = dt.datetime(int(YYYY), int(MM)+1, 1)
    else:
        enddate = dt.datetime(int(YYYY)+1, 1, 1)
    slist = sorted(get_sfile_list(SEISAN_DATA, SEISAN_DB, startdate, enddate))

    for i,sfile in enumerate(slist):
        
        # is this sfile already in summarydf? if so, skip it, unless bool_overwrite
        if summarydf['sfile'].str.contains( os.path.basename(sfile) ).any():
            # already processed
            if not bool_overwrite:
                print('Already processed %s' % sfile )
                continue      
        print('Processing %d of %d: %s' % (i, len(slist), sfile) )
        summarydict = process1event(sfile, bool_overwrite, station_locationsDF=station_locationsDF, MASTER_INV=MASTER_INV)
        summarydf.append(summarydict, ignore_index=True)      
        summarydf.to_csv(summaryfile, index=False)
        filesdone += 1     
        if filesdone >= MAXFILES:
            break  
    return filesdone



def create_details_file(SEISAN_DATA, SEISAN_DB, YYYY, MM, bool_detect_event):
    miniseeddir = os.path.join(SEISAN_DATA, 'miniseed_c', SEISAN_DB)
    summaryfile=os.path.join(miniseeddir, 'summary_%s%s%s.csv' % (SEISAN_DB, YYYY, MM) )
    detailsfile=os.path.join(miniseeddir, 'details_%s%s%s.csv' % (SEISAN_DB, YYYY, MM) )
    if os.path.exists(summaryfile):
        summarydf = pd.read_csv(summaryfile)
        
        for i,row in summarydf.iterrows():
            mseedfile = row['corrected_DSN_mseed']
            if mseedfile:
                corrected_st = read(mseedfile)
                corrected_st.select(component='[ENZ]') # subset to seismic components only
                if len(corrected_st)>0:
                    eventcsvfile = mseedfile.replace('.mseed','.csv')
                    event_df = pd.read_csv(eventcsvfile)
                    eventrow = eventdf2catalogdf(event_df, corrected_st, mseedfile)
                    for key in eventrow.keys():
                        summarydf.loc[i,key] = eventrow[key]
        summarydf.to_csv(detailsfile)          

        
def eventdf2catalogdf(eventdf, st, waveformfilepath, bool_detect_event):
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
    catalogrow['path']=waveformfilepath
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
SEISAN_DATA, SEISAN_DB, station_locationsDF, MASTER_INV, \
bool_detect_event, bool_overwrite = set_globals()
if not isinstance(station_locationsDF, pd.DataFrame):
    print('Station coordinates not loaded. Exiting. Although we could change code and get these from StationXML file')
    exit()
    
# the default is to loop over years/months
filesdone = 0
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
        filesdone = processSeisanYearMonth(SEISAN_DATA, SEISAN_DB, YYYY, MM, filesdone, \
                                           MAXFILES=MAXFILES, bool_overwrite=bool_overwrite, \
                                           station_locationsDF=station_locationsDF, MASTER_INV=MASTER_INV)
        create_details_file(SEISAN_DATA, SEISAN_DB, YYYY, MM, bool_detect_event)

