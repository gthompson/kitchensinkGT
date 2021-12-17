#!/usr/bin/env python
# coding: utf-8
################################
# General setup                #
################################

import os
import sys
import pandas as pd
import pickle
import numpy as np
sys.path.insert(0, './AAA-master/automatic_processing')
#import tools
from config import Config
#from analyzer import Analyzer

import json
from os.path import isfile, isdir
import datetime
from features import FeatureVector
from tools import butter_bandpass_filter
from featuresFunctions import energy, energy_u
from math import sqrt
from tools import print_cm
import time
from tools import extract_features
from copy import deepcopy
import obspy

### PREPARE THE CATALOG DataFrame ###
SEISAN_DATA = os.path.join( os.getenv('HOME'),'DATA','MVO') # e.g. /home/user/seismo
pandaSeisDir = os.path.join(SEISAN_DATA, 'miniseed_c') # e.g. /home/user/seismo/pandaSeis
SEISAN_DB = 'MVOE_' # e.g. the seisan database name (e.g. MVOE_)
PROJECTDIR = os.path.join(os.getenv('HOME'),'src', 'kitchensinkGT', 'PROJECTS', 'MontserratML') # this dir
csvfile_external = os.path.join(SEISAN_DATA, 'MachineLearning', SEISAN_DB, 'runAAA', 'MVOE_11_labelled_events.csv')
csvfile_internal = 'catalog/30_MVO_labelled_events_filtered.csv' # has to match that in AAA-master/config/general/newsettings_10.json
csvfile_internal = './AAA-master/MONTSERRAT/' + csvfile_internal
output_path_cat = csvfile_internal.replace('.csv', '.pd')
alltraces_file = '30_alltraceDFs.csv'

metrics_to_add = ['bandratio_[0.8_4.0_16.0]', 'bandratio_[1.0_6.0_11.0]',                   'bw_max', 'bw_min', 'kurtosis', 'skewness', 'peakF', 'medianF']

# Change if you want your screen to keep quiet
# 0 = quiet
# 1 = in between
# 2 = detailed information
verbatim = 5

# Init project with configuration file
config = Config('./AAA-master/config/general/newsettings_10.json', verbatim=verbatim)
config.readAndCheck()  
cat = pd.read_csv(csvfile_internal)
features = FeatureVector(config, verbatim=verbatim)

# Read all labeled signatures (labels+data) from the catalogue, and extract features
tStart = time.time()
catalog_length = len(cat.index)

# Save featuresList to pickle file
if not os.path.exists('features'):
    os.makedirs('features')
        
# read catalog         
WAVTOPDIR = config.data_to_analyze['path_to_data'] # Glenn. path has miniseed_c hardcoded at start. I want to change this to whatever the config says
for i in range(catalog_length):
    if verbatim > 1:
        print('Processing waveform %d of %d' % (i, catalog_length))
    secondFloat = cat.iloc[i]['second']
    tStartSignature = datetime.datetime(int(cat.iloc[i]['year']),
                                        int(cat.iloc[i]['month']),                                            
                                        int(cat.iloc[i]['day']),                                              
                                        int(cat.iloc[i]['hour']),                                             
                                        int(cat.iloc[i]['minute']),                                           
                                        int(secondFloat),                                         
                                        int((secondFloat-int(secondFloat))*1000000)) #microseconds
    duration = cat.iloc[i]['length']
    path = cat.iloc[i]['path']     
    path = path.replace('miniseed_c', WAVTOPDIR)

    # NEED TO LOAD DATA - FIRST CHECK IF FEATURES ALREADY COMPUTED
    mseedbase = os.path.basename(path)
    print('Trying to read %s' % path)
    if not isfile(path):
        
        # look for WAV file instead, which lacks the .mseed extension
        path_wav = path.replace('.mseed','')
        if not isfile(path_wav): 
            print("File not found: ",path)
            continue
        else:
            st = obspy.read(path_wav)
            fix_trace_id(stream)
    else:
        st = obspy.read(path)
        
    ###### SPIKES CHECK - won't be needed when we reprocess all data #####
    # The reason we do the spike check on the raw WAV file is 
    # because filters run to produce the MSEED file distort the spike
    wavpath = mseedpath.replace('miniseed_c', 'WAV').replace('.mseed','')
    rawst = read(wavpath)
    for tr in rawst:
        check_for_spikes(tr)

    good_traces = 0
    trace_ids_to_eliminate = []
    fix_trace_id(rawst)
    for tr in rawst:
        check_for_spikes(tr)
        if tr.stats.quality_factor > 1.0:
            good_traces += 1
        else:
            trace_ids_to_eliminate.append(tr.id)

    for tr in st:
        if tr.id in trace_ids_to_eliminate:
            st.remove(tr)
    ################ END OF SPIKES CHECK ########################
    
    if len(st)>0 and add_GT_metrics:
        tracecsv = path.replace('.mseed','.csv')
        if isfile(tracecsv):
            tracedf = pd.read_csv(tracecsv)
        else:
            continue
    
    for tr in st:
        featurespkl = os.path.join('features',mseedbase.replace('.mseed', '.%s.pkl' % tr.id)) 
        if os.path.exists(featurespkl):
            continue
        
        # Get signal and its metadata
        s_dict = tr.stats
        
        # Get information about recording
        fs = tr.stats['sampling_rate']         
        length_n = tr.stats['npts'] # only change from read_ubinas
        #d0 = s_dict['starttime']
        #d1 = s_dict['endtime']
        #t_start = datetime.datetime(d0.year,d0.month,d0.day,d0.hour,d0.minute,d0.second)
        #t_end = datetime.datetime(d1.year,d1.month,d1.day,d1.hour,d1.minute,d1.second)
        
        if fs < 70.0 or length_n < 1000:
            continue         

        # Preprocessing & features extraction
        signature = np.array(tr.data)
        print(len(signature), length_n)
        featuresList = extract_features(config, signature, features, fs)
        
        # THIS WOULD BE THE PLACE TO ADD THE PRECOMPUTED FEATURES FROM SEISAN2PANDAS
        if add_GT_metrics:
            thistracedf = tracedf[tracedf['id']==tr.id]
            for col in metrics_to_add:
                featuresList.append(thistracedf[col])
            bandwidth = thistracedf['bw_max'] - thistracedf['bw_min']
            featuresList.append(bandwidth)
        
        with open(featurespkl, 'wb') as f:
            pickle.dump(featuresList, f)

