#!/usr/bin/env python
# coding: utf-8

import os
import shutil
import glob
import pandas as pd
import obspy
from obspy.core import Stream, read 
import numpy as np
from obspy.core.utcdatetime import UTCDateTime
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import datetime
#from sys import exit
from obspy.signal.trigger import classic_sta_lta, z_detect, plot_trigger, trigger_onset
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees
from obspy.taup import TauPyModel

# Glenn Thompson, Feb 2021


#######################################################################
##                Stream tools                                       ##
#######################################################################

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


def removeInstrumentResponse(st, preFilter = (1, 1.5, 30.0, 45.0), outputType = "VEL", inventory = ""):  
    """
    Remove instrument response - assumes inventories have been added to Stream object
    Written for Miami Lakes
    """
    try:
        st.remove_response(output=outputType, pre_filt=preFilter)
    except:
        for tr in st:
            try:
                if inventory:
                    tr.remove_response(output=outputType, pre_filt=preFilter, inventory=inventory)
                else:
                    tr.remove_response(output=outputType, pre_filt=preFilter)
            except:
                print("- Not able to correct data for %s " %  tr.id)
                st.remove(tr)
    return


def eventStatistics(st):
    """ 
    Take a Stream object, measure some basic metrics 
    
    Example: 
        eventStatistics(Stream)

    Returns: 
        list of dictionaries (one per trace)  with keywords id, sample, time, peakamp, energy


    Created for MiamiLakes project.
    
    Note: replace this with method from STA/LTA detection

    """
    # create an empty list of dictionaries
    list_of_dicts = []

    # for each trace, add a new dictionary to the list
    for tr in st:
        thisrow = dict()  # each dictionary will become a row of a dataframe

        # if we use absolute function, we don't need to check if maxy > -miny
        tr.detrend()
        y = np.absolute(tr.data)

        # we did this before
        maxy = y.max()
        maxy_i = y.argmax()
        maxy_t = tr.stats.starttime + maxy_i / tr.stats.sampling_rate

        # add new elements to dictionary
        thisrow['id'] = tr.id
        thisrow['sample'] = maxy_i
        thisrow['time'] = maxy_t
        thisrow['peakamp'] = maxy

        # add a new measurement: energy
        thisrow['energy'] = np.sum(np.square(y)) / tr.stats.sampling_rate

        # add row (dict) to list
        list_of_dicts.append(thisrow)

    # Convert list of dicts to dataframe
    df = pd.DataFrame.from_dict(list_of_dicts)

    return df



def detectEvent(st, pretrig=30, posttrig=30):
    """ 
    Take a Stream object and run an STA/LTA on each channel. Clip out Stream corresponding to first ON trigger and last OFF trigger
    
    Optional inputs: pretrig and posttrig (in seconds).
    
    """
    starttime = st[0].stats.starttime
    mintrig =  999999
    maxtrig = -999999
    for tr in st:
        df = tr.stats.sampling_rate

        print('- detecting')
        cft = z_detect(tr.data, int(10 * df))
        triggerlist = trigger_onset(cft, 0.5, 0.0, max_len = 60 * df)

        for trigpair in triggerlist:
            if trigpair[0] < mintrig:
                mintrig = trigpair[0]
            if trigpair[1] > maxtrig:
                maxtrig = trigpair[1]

    # 30-s before first trigger time or 240 seconds, whichever is greatest as we do not want to trigger more than 1 minute before Steve's times in Excel
    onset_seconds = max([mintrig/df - pretrig, 240])
    # 30-s after last trigger time or 420 seconds, whichever is least as we do not want to trigger more than 2 minutes after Steve's times in Excel
    offset_seconds = min([maxtrig/df + posttrig, 420])
    stime = st[0].stats.starttime + onset_seconds
    etime = st[0].stats.starttime + offset_seconds
    st.trim(starttime = stime, endtime = etime)
    return
    


def mulplt(st, bottomlabel='', ylabels=[]):
    """ Create a plot of a Stream object similar to Seisan's mulplt """
    fh = plt.figure()
    MAXPANELS =6
    n = np.min([MAXPANELS, len(st)])
    
    # start time as a Unix epoch
    startepoch = st[0].stats.starttime.timestamp
    
    # create empty set of subplot handles - without this any change to one affects all
    axh = []
    
    # loop over all stream objects
    for i in range(n):
        # add new axes handle for new subplot
        #axh.append(plt.subplot(n, 1, i+1, sharex=ax))
        axh.append(plt.subplot(n, 1, i+1))
        
        # time vector, t, in seconds since start of record section
        t = np.linspace(st[i].stats.starttime.timestamp - startepoch,
            st[i].stats.endtime.timestamp - startepoch,
            st[i].stats.npts)
            
        # We could detrend, but in case of spikes, subtracting the median may be better
        #st[i].detrend()
        offset = np.median(st[i].data)
        y = st[i].data - offset
        
        # PLOT THE DATA
        axh[i].plot(t, y)


        
    # remove yticks because we will add text showing max and offset values
    axh[i].yaxis.set_ticks([])
        
    # remove xticklabels for all but the bottom subplot
    if i < n-1:
        axh[i].xaxis.set_ticklabels([])
    else:
        # for the bottom subplot, also add an xlabel with start time
        if bottomlabel=='':
            plt.xlabel("Starting at %s" % (st[0].stats.starttime) )
        else:
            plt.xlabel(bottomlabel)
               
    # default ylabel is station.channel
    if ylabels==[]:
        plt.ylabel(st[i].stats.station + "." + st[i].stats.channel, rotation=0)
    else:
        plt.ylabel(ylabels[i])
            
    # explicitly give the maximum amplitude and offset(median)
    plt.text(0, 1, "max=%.1e offset=%.1e" % (np.max(np.abs(y)), offset),
        horizontalalignment='left',
        verticalalignment='top',transform=axh[i].transAxes)
            
    # change all font sizes
    plt.rcParams.update({'font.size': 8})
    
    # show the figure
    plt.show()
    #st.mulplt = types.MethodType(mulplt,st)    
    
    return fh, axh
    
    
def iceweb_spectrogram(st):
    """ Create a plot of a Stream object similar to an IceWeb spectrogram """
    fh = plt.figure()

    # Only do this if want to use a common x-axis
    MAXPANELS =6
    n = np.min([MAXPANELS, len(st)])

    # start time of the record section
    startepoch = st[0].stats.starttime.timestamp

    # create empty set of subplot handles - without this any change to one affects all
    axh = []

    # loop over all stream objects

    for i in range(n):
        # add new axes handle for new subplot
        #axh.append(plt.subplot(n, 1, i+1, sharex=ax))
        axh.append(plt.subplot(2*n, 1, i*2+1))

        # time vector, t, in seconds since start of record section
        t = np.linspace(st[i].stats.starttime.timestamp - startepoch,
                        st[i].stats.endtime.timestamp - startepoch,
                        st[i].stats.npts)

        # We could detrend, but in case of spikes, subtracting the median may be better
        #st[i].detrend()
        offset = np.median(st[i].data)
        y = st[i].data - offset

        # PLOT THE DATA
        axh[i*2+1].plot(t, y)

        # remove yticks because we will add text showing max and offset values
        axh[i*2+1].yaxis.set_ticks([])

        # remove xticklabels for all but the bottom subplot
        if i < n*2-1:
            axh[i].xaxis.set_ticklabels([])
        else:
            # for the bottom subplot, also add an xlabel with wavfilename and start time
            plt.xlabel("WAV file %s\n Starting at %s" % (wavfile, st[0].stats.starttime) )

        # ylabel is station.channel
        plt.ylabel(st[i].stats.station + "." + st[i].stats.channel, rotation=0)

        # explicitly give the maximum amplitude and offset(median)
        plt.text(0, 1, "max=%.1e offset=%.1e" % (np.max(np.abs(y)), offset),
                horizontalalignment='left',
                verticalalignment='top',transform=axh[i].transAxes)

        # PLOT SPECTROGRAM
        axh.append(plt.subplot(n*2, 1, i*2+2))
        #st[i].plot(color='black') # problem here with attaching to subplot
        st[i].spectrogram(axes=axh[i*2+2], log=True, title=st[i].stats.station + " " + str(st[i].stats.starttime))

    # change all font sizes
    plt.rcParams.update({'font.size': 8})

    # show the figure
    plt.show()

    

def max_3c(st):
    """ max of a 3-component seismogram """
    N = len(st)/3
    m = []

    if N.is_integer():
        st.detrend()
        for c in range(int(N)):
            y1 = st[c*3+0].data
            y2 = st[c*3+1].data
            y3 = st[c*3+2].data
            y = np.sqrt(np.square(y1) + np.square(y2) + np.square(y3))
            m.append(max(y))
    return m    
    
    
#######################################################################    
########################         WFDISC tools                        ##
#######################################################################


        
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



#######################################################################
##                BUD tools                                          ##
#######################################################################



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

        this_tr.write(mseedDayFile, format='MSEED') 


    
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
    return


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




######################################################################
##                          FDSN tools                              ##
######################################################################

def get_FDSN_inventory(fdsnClient, eventTime, stationXmlFile, network, latitude, longitude,searchRadiusDeg, pretrigSecs, posttrigSecs ):
    """ 
    Load inventory of stations/channels available around this event time. It will attempt to load from file, then from the client
    Written for Miami Lakes project
    """
    
    if os.path.exists(stationXmlFile):
        # load inv from file
        inv = obspy.core.inventory.read_inventory(stationXmlFile)
    else:
        # load inv from Client & save it        
        startt = eventTime - pretrigSecs
        endt = eventTime + posttrigSecs 
        print('Trying to load inventory from %s to %s' % (startt.strftime('%Y/%m/%d %H:%M'), endt.strftime('%Y/%m/%d %H:%M')))
        try:
            inv = fdsnClient.get_stations(
                network = network,
                latitude = latitude,
                longitude = longitude,
                maxradius = searchRadiusDeg,
                starttime = startt,
                endtime = endt,
                level = 'response'
            )
            
        except Exception as e: 
            print(e)
            print('-  no inventory available')

        else:         
            # Save the inventory to a stationXML file
            inv.write(stationXmlFile, format='STATIONXML')
            print('inventory saved to %s' % stationXmlFile)
            
            return inv



def get_FDSN_Stream(fdsnClient, trace_ids, outfile, startt, endt ):
    """ 
    Load waveform data for all trace ids for this time range
    Written for Miami Lakes project
    """
    if os.path.exists(outfile):
        # load raw data from file
        st = obspy.core.read(outfile)
        
    else:
        # load raw data from FDSN client
        st = obspy.core.Stream()
        for trace_id in trace_ids:
            netw, station, chancode = trace_id.split('.')
            print("net=%s, station=%s, chancode=%s" % (netw, station, chancode))
            try:
                this_st = fdsnClient.get_waveforms(
                    netw,
                    station,
                    "*",
                    chancode,
                    starttime=startt,
                    endtime=endt,
                    attach_response=True
                )

            except:
                print("- No waveform data available for %s for this event %s" % (trace_id, outfile))
                
            else:
                if this_st:
                    this_st.merge(fill_value=0)
                    st += this_st            
    
        if not st:
            print("- No waveform data available for this event %s" % outfile)
            return
        else:
            st.merge(fill_value=0)
    
            # Save raw waveform data to miniseed
            st.write(outfile, format="MSEED")
            
            return st

######################################################################
##                  Inventory tools                                 ##
######################################################################

def inventory2traceid(inv, chancode=''):
    trace_ids = list()

    for networkObject in inv:
        if chancode:
            networkObject = networkObject.select(channel=chancode)
        stationObjects = networkObject.stations

        for stationObject in stationObjects:
            channelObjects = stationObject.channels
            for channelObject in channelObjects:
                this_trace_id = networkObject.code + '.' + stationObject.code + '.' + channelObject.code
                trace_ids.append(this_trace_id)

    return trace_ids


def attach_station_coordinates_from_inventory(inventory, st):
    """ attach_station_coordinates_from_inventory """
    for tr in st:
        for netw in inventory.networks:
            for sta in netw.stations:
                if tr.stats.station == sta.code and netw.code == tr.stats.network:
                    for cha in sta.channels:
                        if tr.stats.location == cha.location_code:
                            tr.stats.latitude = cha.latitude
                            tr.stats.longitude = cha.longitude
    return


######################################################################
##                  Modeling  tools                                 ##
######################################################################


def predict_arrival_times(station, quake):
    """ calculate predicted travel times based on IASP91 model  - see https://docs.obspy.org/packages/obspy.taup.html
        Input: station and quake both are dicts with lat and lon keys
        Output: a phases dict is added to station, with phase name keys and predicted arrival times """
    model = TauPyModel(model="iasp91")
    
    [dist_in_m, az1, az2] = gps2dist_azimuth(quake['lat'], quake['lon'], station['lat'], station['lon'])
    station['distance'] = kilometers2degrees(dist_in_m/1000)
    arrivals = model.get_travel_times(source_depth_in_km=quake['depth'],distance_in_degree=station['distance'])
    # https://docs.obspy.org/packages/autogen/obspy.taup.helper_classes.Arrival.html#obspy.taup.helper_classes.Arrival
    
    phases = dict()
    for a in arrivals:
        phasetime = quake['otime'] + a.time
        phases[a.name] = phasetime.strftime('%H:%M:%S')
        if a.name == 'S':
            Rtime = quake['otime'] + a.time/ ((0.8453)**0.5)
            phases['Rayleigh'] = Rtime.strftime('%H:%M:%S')
    station['phases'] = phases
    
    return station

def syngine2stream(station, lat, lon, GCMTeventID, mseedfile):
    """ Generate synthetics for a GCMT event, save into an mseedfile, return as Stream object """
    if os.path.exists(mseedfile):
        synth_disp = read(mseedfile)
    else:
        synth_disp = read("http://service.iris.edu/irisws/syngine/1/query?"
                  "format=miniseed&units=displacement&dt=0.02&"
                  "receiverlatitude=%f&receiverlongitude=%f&"
                  "eventid=GCMT:%s" % (lat, lon, GCMTeventID))
        for c in range(len(synth_disp)):
            synth_disp[c].stats.latitude = lat
            synth_disp[c].stats.longitude = lon
        synth_disp.write(mseedfile)
    return synth_disp



