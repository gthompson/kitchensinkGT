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
                           bool_overwrite=False, station_locationsDF=None, MASTER_INV=None, bool_index_only=False):
    
    
    # The summary file created here is just an index to different files
    miniseeddir = os.path.join(SEISAN_DATA, 'miniseed_c', SEISAN_DB)
    sfileindexfile=os.path.join(miniseeddir, 'sfileindex_%s%s%s.csv' % (SEISAN_DB, YYYY, MM) )

    if os.path.exists(sfileindexfile) and bool_overwrite:
        sfileindex_df = pd.read_csv(summaryfile)
    else:
        sfileindex_df = pd.DataFrame(columns=['sfile', 'DSN_wavfile', 'DSN_exists', 'ASN_wavfile', 'ASN_exists', 'corrected_DSN_mseed', 'corrected_ASN_mseed', 'mainclass', 'subclass']) 

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
        if sfileindex_df['sfile'].str.contains( os.path.basename(sfile) ).any():
            # already processed
            if not bool_overwrite:
                print('Already processed %s' % sfile )
                continue      
        print('Processing %d of %d: %s' % (i, len(slist), sfile) )
        try:
            sfileindex_dict = process1event(sfile, bool_overwrite, station_locationsDF=station_locationsDF, \
                                    MASTER_INV=MASTER_INV, bool_index_only=bool_index_only)
            sfileindex_df.append(sfileindex_dict, ignore_index=True)  
        except:
            if len(sfileindex_df.index)>1:
                sfileindex_df.to_csv(sfileindexfile, index=False)
        filesdone += 1     
        if filesdone >= MAXFILES:
            break  
    sfileindex_df.to_csv(sfileindexfile, index=False)        
    return filesdone



#######################################################################
print('Usage %s STARTYEAR STARTMONTH MAXFILES bool_index_only bool_overwrite' % sys.argv[0])
    
# import controlling variables - actually not globals
SEISAN_DATA, SEISAN_DB, station_locationsDF, MASTER_INV, bool_overwrite = set_globals()
if not isinstance(station_locationsDF, pd.DataFrame):
    print('Station coordinates not loaded. Exiting. Although we could change code and get these from StationXML file')
    exit()    
STARTYEAR = "1996"
STARTMONTH = "01"
MAXFILES = 999999
bool_index_only = False

# command line arguments
argc = len(sys.argv)
print(argc)
if argc>1:
    STARTYEAR = sys.argv[1]
if argc>2:
    STARTMONTH = sys.argv[2]
if argc>3:
    MAXFILES = int(sys.argv[3])
if argc>4:
    if sys.argv[4]=='1':
        bool_index_only = True
if argc>5:
    if sys.argv[5]=='1':
        bool_overwrite = True  

    
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
                                           station_locationsDF=station_locationsDF, \
                                           MASTER_INV=MASTER_INV, bool_index_only=bool_index_only)


