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
from libMVO import fix_trace_id, inventory_fix_id_mvo, load_mvo_inventory
from metrics import process_trace, choose_best_traces, select_by_index_list, ampengfft
from libseisGT import Stream_min_starttime, detect_network_event
from seisan_classes import spath2datetime, Sfile #, printEvents

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




def wav2stream(paths):
    if os.path.exists(paths['picklefile_c']):
        st = read(paths['picklefile_c'])
        return st
    if os.path.exists(paths['picklefile_u']):
        st = read(paths['picklefile_u'])
    else:
        # We need to create the pickle file
        # Try to read the WAV file
        st = Stream()
        print('Processing %s.' % paths['wavbase'], end=' ')
        print('Reading.', end = ' ')
        try:               
            st = read(paths['wavfile'])
        except:
            print('ERROR. Could not load.')
            return st

        if len(st)==0:
            print('ERROR. No traces.')
            return st
        else: 
            print('Success.')
        
        fix_trace_id(st, shortperiod=bool_shortperiod) 
        #st=st.select(component='Z')

        for tr in st:
            this_inv = None            
            process_trace(tr, inv=this_inv)
        
        # remove bad traces
        for tr in st:    
            if tr.stats.quality_factor <= 0.0:
                st.remove(tr)  
                
            
        # Write pickle file. This is now the best place to load stream data from in future, so replaces
        # seisan/WAV directory with seisan/PICKLE
        if len(st)>0:
            print('Writing ',paths['picklefile_u'])
            st.write(paths['picklefile_u'], format='PICKLE') 

    for tr in st:
        this_inv = None            
        if bool_correct_data: # try to find corresponding station XML
            this_inv = load_mvo_inventory(tr, paths['CALDIR'])
        process_trace(tr, inv=this_inv)
    if len(st)>0:
        print('Writing ',paths['picklefile_c'])
        st.write(paths['picklefile_c'], format='PICKLE') 
        
    return st

def Stream2logfile(st, paths):        
    # save log file
    if not os.path.exists(paths['logfile']):
        print('Writing %s' % paths['logfile'])
        with open(paths['logfile'],'w') as fout:
            for tr in st:
                # print trace history
                pprint(tr.stats, stream=fout)
                fout.write('\n')
                
def Stream2png(st, paths):

    if not os.path.exists(paths['seismicpngfile']):
        print('Writing ',paths['seismicpngfile'])                
        chosen = choose_best_traces(st, MAX_TRACES=99, include_seismic=True, 
                                    include_infrasound=False, include_uncorrected=False)
        if len(chosen)>0:
            st2 = select_by_index_list(st, chosen)
            st2.plot(equal_scale=False, outfile=paths['seismicpngfile'])   
            
    if not os.path.exists(paths['infrasoundpngfile']):
        print('Writing ',paths['infrasoundpngfile'])                
        chosen = choose_best_traces(st, MAX_TRACES=99, include_seismic=False, 
                                    include_infrasound=True, include_uncorrected=False)
        if len(chosen)>0:
            st2 = select_by_index_list(st, chosen)
            st2.plot(equal_scale=False, outfile=paths['infrasoundpngfile'])    
            
def add_ampengfft_metrics(st, paths):
    print('Computing spectrogram data')        
    iwsobj = IceWeb.icewebSpectrogram(stream=st)
    iwsobj = iwsobj.precompute() # spectrograms data added
    print('Computing spectrum.', end = ' ')  
    iwsobj.compute_amplitude_spectrum(compute_bandwidth=True) # adds tr.stats.spectrum
    for tr in iwsobj.stream:
        ampengfft(tr, paths['PICKLEDIR']) # add peaktime, peakamp, energy
    return iwsobj


def plot_spectrograms(paths, iwsobj):
    st = iwsobj.stream

    # free scale  
    print('Creating %s.' % paths['sgramfile'], end = ' ')            
    titlestr = os.path.basename(paths['sgramfile']) 
    chosen = choose_best_traces(st, MAX_TRACES=10)            
    iwsobj.plot(outfile=paths['sgramfile'], log=False, equal_scale=False, add_colorbar=True, dbscale=True, title=titlestr, trace_indexes=chosen);

    # fixed scale 
    print('Creating %s.' % paths['sgramfixed'], end = ' ')
    titlestr = os.path.basename(paths['sgramfixed'])
    clim_in_dB = [-160, -100] # only works for corrected volcano-seismic data
    clim_in_units = [ IceWeb.dB2amp(clim_in_dB[0]),  IceWeb.dB2amp(clim_in_dB[1]) ]
    iwsobj.plot(outfile=paths['sgramfixed'], log=False, clim=clim_in_units, add_colorbar=True, dbscale=True, title=titlestr, trace_indexes=chosen);           
    # need separate plots for any infrasound channels or uncorrected data

def metrics2tracedf(traceCSVfile, st):
    if os.path.exists(traceCSVfile):
        tracedf = pd.read_csv(traceCSVfile)
    else: 
        print('- Building metrics dataframe.', end = ' ') 
        tracedf = pd.DataFrame()
        list_of_tracerows = []
        for tr in st:
            s = tr.stats
            tracerow = {'id':tr.id, 'starttime':s.starttime, 
                   'Fs':s.sampling_rate, 
                   'calib':s.calib, 'units':s.units, 
                   'quality':s.quality_factor}
            if 'spectrum' in s: 
                for item in ['medianF', 'peakF', 'peakA', 'bw_min', 'bw_max']:
                    try:
                        tracerow[item] = s.spectrum[item]
                    except:
                        pass
            if 'metrics' in s:
                m = s.metrics
                for item in ['snr', 'signal_level', 'noise_level', 'twin',
                             'peakamp', 'peaktime', 'energy', 'RSAM_high', 'RSAM_low',
                             'sample_min', 'sample_max', 'sample_mean', 'sample_median', 
                             'sample_lower_quartile', 'sample_upper_quartile', 'sample_rms', 
                             'sample_stdev', 'percent_availability', 'num_gaps', 'skewness', 'kurtosis']:
                             #'start_gap', 'num_gaps', 'end_gap', 'sum_gaps', 'max_gap', 
                             #'num_overlaps', 'sum_overlaps', 'num_records', 'record_length', 
                    try:
                        tracerow[item] = m[item]
                    except:
                        pass 
            if 'bandratio' in s:
                for dictitem in s['bandratio']:
                    label = 'bandratio_' +  "".join(str(dictitem['freqlims'])).replace(', ','_')
                    tracerow[label] = dictitem['RSAM_ratio']

            list_of_tracerows.append(tracerow)
        tracedf = pd.DataFrame(list_of_tracerows)
        tracedf = tracedf.round({'Fs': 2, 'secs': 2, 'quality':2, 'medianF':1, 'peakF':1, 'bw_max':1, 'bw_min':1, 'peaktime':2, 'twin':2, 'skewness':2, 'kurtosis':2})
        print('Saving to CSV.')
        tracedf.set_index('id')
        tracedf.to_csv(traceCSVfile)
    return tracedf

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
    paths['HTMLDIR'] = paths['WAVDIR'].replace('WAV', 'HTML')
    paths['PICKLEDIR'] = paths['WAVDIR'].replace('WAV', 'PICKLE')           
    paths['seismicpngfile'] = os.path.join(paths['HTMLDIR'], paths['wavbase'] + '_seismic.png')
    paths['infrasoundpngfile'] = os.path.join(paths['HTMLDIR'], paths['wavbase'] + '_infrasound.png')
    paths['picklefile_u'] = os.path.join(paths['PICKLEDIR'], paths['wavbase'] + '_u.pickle')
    paths['picklefile_c'] = os.path.join(paths['PICKLEDIR'], paths['wavbase'] + '_c.pickle')
    paths['logfile'] = paths['picklefile'].replace('.pickle', '.log')
    paths['sgramfile'] = os.path.join(paths['HTMLDIR'], paths['wavbase'] + '_sgram.png')
    paths['sgramfixed'] = paths['sgramfile'].replace('_sgram.png', '_sgram_fixed.png')  
    paths['traceCSVfile_u'] = paths['picklefile_u'].replace('.pickle', '.csv')
    paths['traceCSVfile_c'] = paths['picklefile_c'].replace('.pickle', '.csv')
    paths['detectionfile'] = os.path.join(paths['HTMLDIR'], paths['wavbase'] + '_detection.png')

    return paths
    
    
def examineWAV(wavfile):
    paths = wavfile2paths(wavfile)
    if not os.path.exists(paths['HTMLDIR']) and bool_make_png_files:
        os.makedirs(paths['HTMLDIR'])
    if not os.path.exists(paths['PICKLEDIR']):
        os.makedirs(paths['PICKLEDIR'])  
    if bool_overwrite:
        products = ['picklefile_u', 'picklefile_c', 'seismicpngfile', 'infrasoundpngfile', 'logfile', 'sgramfile', 
                    'sgramfixed', 'traceCSVfile_u', 'traceCSVfile_c', 'detectionfile']
        for product in products:
            if os.path.exists(paths[product]):
                os.remove(paths[product])
    wavrow={}
    
    st = wav2stream(paths)
    st.select(component='[ENZ]') # subset to seismic components only
    if len(st)==0:
        return wavrow
    Stream2logfile(st, paths)
    if bool_make_png_files:
        Stream2png(st, paths)
        
    if not 'energy' in st[0].stats.metrics: 
        iwsobj = add_ampengfft_metrics(st, paths)
        if not os.path.exists(paths['sgramfile']) and bool_make_png_files:
            plot_spectrograms(paths, iwsobj)
        for tr in iwsobj.stream:
            tr.stats.pop('spectrogramdata', None) # Remove spectrogramdata as it is too large for picklefile
            #tr.stats.spectrogramdata = {}
        print('Writing enhanced pickle file.')     
        st = iwsobj.stream
        st.write(paths['picklefile'], 'PICKLE') # Rewrite pickle file with extra attributes     
        
    tracedf = metrics2tracedf(paths['traceCSVfile'], st)
    wavrow = tracedf2eventdf(tracedf, st, paths)        
    return wavrow


def examineSeisanYearMonth(SEISAN_DATA, DB, YYYY, MM, filesdone, MAXFILES=999999, startdate=None):
    
    # Get s-file list
    if not startdate:
        startdate = dt.datetime(int(YYYY), int(MM), 1)
    if int(MM)<12:
        enddate = dt.datetime(int(YYYY), int(MM)+1, 1)
    else:
        enddate = dt.datetime(int(YYYY)+1, 1, 1)
    slist = sorted(get_sfile_list(SEISAN_DATA, DB, startdate, enddate))
    
    for i,sfile in enumerate(slist):
        print('Processing %d of %d: %s' % (i, len(slist), sfile) )
        if i==MAXFILES:
            break
            
        s = Sfile(sfile, use_mvo_parser=True)
        #s.cat()
        #s.printEvents()
        d = s.to_dict()
        #pprint(d)
   	
        # SCAFFOLD - Add code here to check wavfiles are consistent with S-file name. I have this in a later Jupyter Notebook.
     
        for item in ['wavfile1', 'wavfile2']:
            if d[item]:
                if os.path.exists(d[item]):
                    wavbase = os.path.basename(d[item])
                    if 'MVO' in wavbase:
                        print('Processing ',d[item])
                        examineWAV(d[item])

def main(DB, MAXFILES=999999, startdate=None):
    filesdone = 0
    yeardirs = sorted(glob(os.path.join('REA',DB,'[12]???')))
    for yeardir in yeardirs:
        YYYY = os.path.basename(yeardir)
        monthsdirs = sorted(glob(os.path.join(yeardir,'[01]?')))
        for monthdir in monthsdirs:
            if filesdone>=MAXFILES:
                break
            MM = os.path.basename(monthdir)
            print('**** Processing %s ****' % monthdir)
            examineSeisanYearMonth('.', DB, YYYY, MM, filesdone, MAXFILES=999999, startdate=startdate)

SEISAN_DATA = os.path.join( os.getenv('HOME'),'DATA','MVO')
os.chdir(SEISAN_DATA)
SEISAN_DB = 'MVOE_'
bool_shortperiod=False
bool_correct_data=True
bool_make_png_files=False
bool_detect_event=True
bool_overwrite=False
startdate = dt.datetime(1995,7,18,0,0,0)
startdate = dt.datetime(2000,3,1,0,0,0)
# 2000-02-29-2338-36S.MVO___019
main(SEISAN_DB, 999999, startdate)
