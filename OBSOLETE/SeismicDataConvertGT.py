#!/usr/bin/env python
# coding: utf-8


# Library of Python functions created to help convert ROCKETSEIS and CALIPSO data


import os
print(os.getcwd())
import shutil
import glob
import pandas as pd
from obspy.core import Stream, read 
import numpy as np
from obspy.core.utcdatetime import UTCDateTime
import warnings
warnings.filterwarnings("ignore")



def decompressFile(compressedFilePath, decompressDir):
    """ 
    Gunzip and/or untar a compressed file.
    
    Many of the data files are in tar.Z format.


    Created for CALIPSO data archive from Alan Linde.   
    """

    pwd = os.getcwd()
    compressedFilePath = os.path.abspath(compressedFilePath)
    decompressDir = os.path.abspath(decompressDir)
    
    bname = os.path.basename(compressedFilePath)
    if os.path.exists(compressedFilePath):
        if not os.path.exists(decompressDir):
            os.makedirs(decompressDir)        
        newcompressedFilePath = os.path.join(decompressDir, bname)
        shutil.copy(compressedFilePath, newcompressedFilePath)
        os.chdir(decompressDir)
    else:
        print('No such file: %s' % compressedFilePath)
        return
    
    # Check for Zip file
    parts = bname.split('.')
    if parts[-1].lower().find('z') > -1:
        
        # zipfile
        print('Trying gunzip: %s' % bname)
        try:
            os.system('gunzip %s' % bname)
            print('- success')
            os.system('mv %s %s.PROCESSED' % (compressedFilePath, compressedFilePath))
        except:
            print('not a gzip file %s' % bname)
            os.chdir(pwd)
            return
        del parts[-1]
        sep = '.'
        bname = sep.join(parts)
    else:
        print('Not a zipfile %s' % bname)
        
    # Check for TAR file
    if parts[-1].lower().find('tar') > -1:
        if os.path.exists(bname):
            print('Trying to untar %s' % bname)
            try:
                os.system('tar -xf %s' % bname)
                print('- success')
                os.system('rm %s' % bname)
                if os.path.exists(compressedFilePath):
                    os.system('mv %s %s.PROCESSED' % (compressedFilePath, compressedFilePath))
            except:
                print('Cannot untar %s' % bname)
        else:
            print('Not found %s' % bname)
    else:
        print('Not a tar file %s' % bname)
    
    
    os.chdir(pwd)
    return





def Stream_min_starttime(all_traces):
    """
    Take a Stream object, and return the minimum starttime

    Created for CALIPSO data archive from Alan Linde.
    """ 

    min_stime = UTCDateTime(2099, 12, 31, 0, 0, 0.0)
    max_stime = UTCDateTime(1900, 1, 1, 0, 0, 0.0)
    min_etime = UTCDateTime(2099, 12, 31, 0, 0, 0.0)
    max_etime = UTCDateTime(1900, 1, 1, 0, 0, 0.0)    
    for this_tr in all_traces:
        if this_tr.stats.starttime < min_stime:
            min_stime = this_tr.stats.starttime
        if this_tr.stats.starttime > max_stime:
            max_stime = this_tr.stats.starttime  
        if this_tr.stats.endtime < min_etime:
            min_etime = this_tr.stats.endtime
        if this_tr.stats.endtime > max_etime:
            max_etime = this_tr.stats.endtime              
    return min_stime, max_stime, min_etime, max_etime



def Stream_to_24H(all_traces):
    """
    Take a Stream object, merge all traces with common ids and pad out to 24-hour-long traces

    Created for ROCKETSEIS data conversion and modified for CALIPSO data archive from Alan Linde. 
    """

    all_traces.merge(fill_value=0)
    min_stime, max_stime, min_etime, max_etime = Stream_min_starttime(all_traces)
    
    desired_stime = UTCDateTime(min_stime.year, min_stime.month, min_stime.day, 0, 0, 0.0)
    desired_etime = desired_stime + 86400
    
    days = Stream()
    while True:
        
        this_st = all_traces.copy()
        this_st.trim(starttime=desired_stime, endtime=desired_etime, pad=True, fill_value=0)
        for this_tr in this_st:
            days.append(this_tr)
        desired_stime += 86400
        desired_etime += 86400
        if desired_etime > max_etime + 86400:
            break
    return days



def Trace_merge_with_BUDfile(this_tr, budfile):
    """
    Clever way to merge overlapping traces into a BUD file. Uses all non-zero data values from both.

    Created for CALIPSO data archive from Alan Linde, when needed to upgrade Stream_to_BUD.
    """

    other_st = read(budfile)
    error_flag = False
    
    if len(other_st)>1:
        print('More than 1 trace in %s. Cannot merge.' % budfile)
        error_flag = True
        
    other_tr = other_st[0]
    if not (this_tr.id == other_tr.id):
        print('Different trace IDs. Cannot merge.')
        error_flag = True
        
    if not (this_tr.stats.sampling_rate == other_tr.stats.sampling_rate):
        print('Different sampling rates. Cannot merge.')
        error_flag = True
        
    if (abs(this_tr.stats.starttime - other_tr.stats.starttime) > this_tr.stats.delta/4):
        print('Different start times. Cannot merge.')  
        error_flag = True

    if (abs(this_tr.stats.endtime - other_tr.stats.endtime) > this_tr.stats.delta/4):
        print('Different end times. Cannot merge.')  
        error_flag = True
        
    if error_flag: # traces incompatible, so return the trace with the most non-zero values
        this_good = np.count_nonzero(this_tr.data)
        #print(this_tr.stats)
        other_good = np.count_nonzero(other_tr.data)
        #print(other_tr.stats)
        if other_good > this_good:
            return other_tr
        else:
            return this_tr
    
    else: # things are good
        indices = np.where(other_tr.data == 0)
        other_tr.data[indices] = this_tr.data[indices]
        return other_tr
    



def Stream_to_BUD(TOPDIR, all_traces):
    """ 
    Take a Stream object and write it out in IRIS/PASSCAL BUD format. 
    
    Example:
    
        Stream_to_BUD('RAW', all_traces)
    Where all_traces is a Stream object with traces from 2020346

    Creates a BUD directory structure that looks like:

        DAYS
        ├── BHP2
        │   ├── 1R.BHP2..EH1.2020.346
        │   ├── 1R.BHP2..EH2.2020.346
        │   └── 1R.BHP2..EHZ.2020.346
        ├── BHP4
        │   ├── 1R.BHP4..EH1.2020.346
        │   ├── 1R.BHP4..EH2.2020.346
        │   └── 1R.BHP4..EHZ.2020.346
        ├── FIREP
        │   ├── 1R.FIREP..EH1.2020.346
        │   ├── 1R.FIREP..EH2.2020.346
        │   └── 1R.FIREP..EHZ.2020.346
        └── TANKP
            ├── 1R.TANKP..EH1.2020.346
            ├── 1R.TANKP..EH2.2020.346
            └── 1R.TANKP..EHZ.2020.346
        
    where BHP2, BHP4, FIREP and TANKP are station names, 1R is network name, 
    location is blank, channels are EH[Z12], year is 2020 and day of year is 346.

    Created for ROCKETSEIS data conversion and modified for CALIPSO data archive from Alan Linde.   
    """
    
    all_traces = Stream_to_24H(all_traces)
    
    daysDir = os.path.join(TOPDIR, 'DAYS')

    for this_tr in all_traces:
        YYYY = this_tr.stats.starttime.year
        JJJ = this_tr.stats.starttime.julday
        stationDaysDir = os.path.join(daysDir, this_tr.stats.station)
        if not os.path.exists(stationDaysDir):
            os.makedirs(stationDaysDir)
            #print(stationDaysDir)
        mseedDayBasename = "%s.%04d.%03d" % (this_tr.id, YYYY, JJJ  )
        mseedDayFile = os.path.join(stationDaysDir, mseedDayBasename)
        #print(mseedDayFile)
        if os.path.exists(mseedDayFile):
            this_tr = Trace_merge_with_BUDfile(this_tr, mseedDayFile)
        print('Writing ',this_tr,' to ',mseedDayFile)
        this_tr.write(mseedDayFile, format='MSEED') 





        
def index_waveformfiles(wffiles):
    """ 
    Take a list of seismic waveform data files and return a dataframe similar to a wfdisc table

    Created for CALIPSO data archive from Alan Linde.
    """

    wfdisc_df = pd.DataFrame()
    traceids = []
    starttimes = []
    endtimes = []
    sampling_rates = []
    calibs = []
    ddirs = []
    dfiles = []
    npts = []
    for wffile in sorted(wffiles):
        dfile = os.path.basename(wffile)
        ddir = os.path.dirname(wffile)
        try:
            this_st = read(wffile)
            print('Read %s\n' % wffile)
        except:
            print('Could not read %s\n' % wffile)
            next
        else:
            for this_tr in this_st:
                r = this_tr.stats
                traceids.append(this_tr.id)
                starttimes.append(r.starttime)
                endtimes.append(r.endtime)
                sampling_rates.append(r.sampling_rate)
                calibs.append(r.calib)
                ddirs.append(ddir)
                dfiles.append(dfile)
                npts.append(r.npts)
    if wffiles:
        wfdisc_dict = {'traceID':traceids, 'starttime':starttimes, 'endtime':endtimes, 'npts':npts, 
                       'sampling_rate':sampling_rates, 'calib':calibs, 'ddir':ddirs, 'dfile':dfiles}
        #print(wfdisc_dict)
        wfdisc_df = pd.DataFrame.from_dict(wfdisc_dict)  
        wfdisc_df.sort_values(['starttime'], ascending=[True], inplace=True)
    return wfdisc_df



def wfdisc_to_BUD(wfdisc_df, TOPDIR, put_away):
    """ 
    Read a wfdisc-like dataframe, read the associated files, and write out as BUD format


    Created for CALIPSO data archive from Alan Linde.
    """

    unique_traceIDs = wfdisc_df['traceID'].unique().tolist()
    print(unique_traceIDs)
    
    successful_wffiles = list()

    for traceID in unique_traceIDs:
        print(traceID)
        
        trace_df = wfdisc_df[wfdisc_df['traceID']==traceID]
        
        # identify earliest start time and latest end time for this channel
        #print(trace_df.iloc[0]['starttime'])
        #print(trace_df.iloc[-1]['endtime'])
        minUTC = trace_df.starttime.min()
        maxUTC = trace_df.endtime.max()
        start_date = minUTC.replace(hour=0, minute=0, second=0, microsecond=0)
        end_date = maxUTC.replace(hour=23, minute=59, second=59, microsecond=999999)
        this_date = start_date

        while this_date <= end_date: 
            all_traces = Stream()
        
            # loop from earliest start day to latest end day
            subset_df = trace_df[(trace_df['starttime'] < this_date+86400) & (trace_df['endtime'] >= this_date)]
            #print(subset_df)
            
            if len(subset_df.index)==0:
                next
        
            for index, row in subset_df.iterrows():
                wffile = os.path.join(row['ddir'], row['dfile'])
                start_at = max([this_date, row['starttime']])
                end_at = min([this_date+86400, row['endtime']])
                print('- ',wffile,': START AT:', start_at, ', END AT: ',end_at)
                try:
                    this_st = read(wffile, starttime=start_at, endtime=end_at)

                except:
                    print(' Failed\n')
                    next
                else:
                    print(' Succeeded\n')
                    #raise Exception("Stopping here")
                    if end_at == row['endtime']:
                        successful_wffiles.append(wffile)   
                    for this_tr in this_st:
                        if this_tr.id == traceID:
                            #print(tr.stats)
                            all_traces = all_traces.append(this_tr)
            #print(st.__str__(extended=True)) 
            try:
                all_traces.merge(fill_value=0)
            except:
                print('Failed to merge ', all_traces)
            print(all_traces.__str__(extended=True))
        
            # Check that we really only have a single trace ID before writing the BUD files
            error_flag = False
            for this_tr in all_traces:
                if not this_tr.id == traceID:
                    error_flag = True
            if not error_flag:
                try:
                    Stream_to_BUD(TOPDIR, all_traces)
                except:
                    print('Stream_to_BUD failed for ', all_traces)
            
            this_date += 86400
            
    for wffile in successful_wffiles:
        ddir = os.path.dirname(wffile)
        dbase = "%s.PROCESSED" % os.path.basename(wffile)
        newwffile = os.path.join(ddir, dbase)
        print('move %s %s' % (wffile, newwffile))
        if os.path.exists(wffile) and put_away:
            shutil.move(wffile, newwffile)

            
            
def process_wfdirs(wfdirs, filematch, put_away=False):
    """ 
    Process a directory containing waveform data files in any format readable by ObsPy.
    Build a wfdisc-like dataframe indexing those waveform files.
    Convert them to a BUD archive.

    Created for CALIPSO data archive from Alan Linde.
    """

    for wfdir in wfdirs:
        print('Processing %s' % wfdir)
        wffiles = glob.glob(os.path.join(wfdir, filematch))
        if wffiles:
            #print(wffiles)
            wfdisc_df = index_waveformfiles(wffiles)
            #print(wfdisc_df)
            if not wfdisc_df.empty:
                wfdisc_to_BUD(wfdisc_df, TOPDIR, put_away)  
    print('Done.')


    
def BUD_load_day(BUDDIR, year, jday):
    """
    Load all files corresponding to this year and day from a BUD archive


    Created for CALIPSO data archive from Alan Linde.
    """

    all_stations = glob.glob(os.path.join(BUDDIR, '*'))
    all_traces = Stream()
    for station_dir in all_stations:
        all_files = glob.glob(os.path.join(station_dir, '*.%04d.%03d' % (year, jday)))
        for this_file in all_files:
            try:
                these_traces = read(this_file)
            except:
                print('Cannot read %s' % this_file)
            else:
                for this_tr in these_traces:
                    all_traces.append(this_tr)
    return all_traces


def rt130YYYYJJJdirectory_to_StreamDay(TOPDIR, YYYYJJJ, NETWORK, chans, digitizerStationPairs):
    """ 
    Explore a YYYYJJJ directory from one or many RT130 digitizers. 
    Read the Reftek waveform data files into an ObsPy Stream object that is 24-hours long.
    
    A YYYYJJJ directory has a directory structure like:
    
        2020346
        ├── 91F8
        │   ├── 0
        │   ├── 1
        │   └── 9
        ├── 9338
        │   ├── 0
        │   └── 1
        ├── 9D7C
        │   ├── 0
        │   └── 1
        └── AB13
            ├── 0
            ├── 1
            └── 9
    
    where 91F8, 9338, 9D7C and AB13 are the names of Reftek RT130 digitizers, 
    and the 1 directories contain Reftek waveform data in (usually) 1-hour chunks.
    
    Inputs:
        TOPDIR                - the directory beneath which the YYYYJJJ directories are stored
        YYYYJJJ               - the YYYYJJJ directory to process
        NETWORK               - the 2-character seismic network ID
        chans                 - a list of channel names for each station, usually EHZ, EH1 and EH2
        digitizerStationPairs - some Reftek waveform files do not contain station metadata in the header
                                Missing pairings can be described in a dictionary

    Example:
        TOPDIR = './RAW'
        YYYYJJJ = '2020346'
        NETWORK = '1R'
        chans = ['EHZ','EH1','EH2']
        digitizerStationPairs = {'9338':'BHP4'}    
        all_traces = rt130YYYYJJJdirectory_to_StreamDay(TOPDIR, YYYYJJJ, NETWORK, chans, digitizerStationPairs)


    Created to convert ROCKETSEIS data.    
    """

    all_traces = Stream()
    digitizerDirs = glob.glob(os.path.join(TOPDIR, YYYYJJJ, '????'))

    # Convert the 1 directories
    for thisDigitizerDir in digitizerDirs:
        oneDir = os.path.join(thisDigitizerDir, '1')
        if os.path.exists(oneDir):
            rt130Files = sorted(glob.glob(os.path.join(oneDir, '?????????_????????')))
            for rt130File in rt130Files:
                try:
                    st = read(rt130File)  
                    for c in range(3):
                        st[c].stats.network = NETWORK
                        st[c].stats.channel = chans[c]
                except:
                    print(rt130File, ' is bad\n')                
                else:

                    for tr in st:
                        if not tr.stats.station:
                            print('No station metadata in %s' % rt130File)
                            thisDigitizer = os.path.basename(thisDigitizerDir)
                            if digitizerStationPairs[thisDigitizer]:
                                tr.stats.station = digitizerStationPairs[thisDigitizer]

                        if tr.stats.station:
                            try:
                                all_traces.append(tr)
                            except:
                                print("Could not add trace")
                                
    all_traces = Stream_to_24H(all_traces)
    
    return all_traces


def Stream_to_dayplot(TOPDIR, all_traces):
    """ 
    Take a Stream object, pad it to 24-hours, plot it, and save to a PNG file. 
    
    Example: 
        Stream_to_dayplot('RAW', all_traces)

    Creates: 
    
        DAYS
        ├── 1R.2020.346.png


    Make sure that all_traces[0] contains full trace-id metadata. 

    Created for ROCKETSEIS project.

    """    

    daysDir = os.path.join(TOPDIR, 'DAYPLOTS')
    os.makedirs(daysDir)
    NETWORK = all_traces[0].stats.network
    stime = all_traces[0].stats.starttime
    YYYY = stime.year
    JJJ = stime.yearday
    pngfile = os.path.join(daysDir, '%s.%s.%s.png' % (NETWORK, YYYYJJJ[0:4], YYYYJJJ[4:])  )   
    all_traces.plot(equal_scale=False, outfile=pngfile);

