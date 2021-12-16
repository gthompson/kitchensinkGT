# -*-coding:Utf-8 -*

# Copyright: Marielle MALFANTE - GIPSA-Lab -
# Univ. Grenoble Alpes, CNRS, Grenoble INP, GIPSA-lab, 38000 Grenoble, France
# (04/2018)
#
# marielle.malfante@gipsa-lab.fr (@gmail.com)
#
# This software is a computer program whose purpose is to automatically
# processing time series (automatic classification, detection). The architecture
# is based on machine learning tools.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL

import obspy
import soundfile
import numpy as np
from os.path import isfile
import matplotlib.pylab as plt
import datetime
from obspy import UTCDateTime
from obspy.clients.arclink.client import Client
from tools import filter_data

# added by Glenn
import pandas as pd 
import os, sys
LIBpath = os.path.join( os.getenv('HOME'),'src','kitchensinkGT', 'LIB')
sys.path.append(LIBpath)
from libMVO import fix_trace_id

def read_ubinas(file_path, config, verbatim=0):
    """ Read Ubinas data.
    /!\ Signature of this function should not be modified and is similar for all applications (i,e, read_XXX)
    INPUT:
    - file_path: to file containing data (.dat here)
    - config: config dictionnary according to project formating
    - verbatim
    OUTPUT:
    - s: numpy array containing the signal read
    - fs: sampling_rate
    - t_start and t_end: as datetime.datetime objects
    - length_n
    """

    # Read and check reading
    if not isfile(file_path):
        print("No file at",file_path)
        return None, None
    stream = obspy.read(file_path)
    if len(stream.traces) == 0:
        if verbatim > 0:
            print('Could not read ',file_path)
        return None
    trace = stream.traces[0]
    # Get signal and its metadata
    s = trace.data
    s_dict = trace.stats
    # Gets information about recording
    fs = s_dict['sampling_rate']
    length_n=s_dict['sac']['npts']
    if length_n != np.shape(s)[0]:
        print('problem with signals reading dimensions:', length_n, np.shape(s)[0])
    d0 = s_dict['starttime']
    d1 = s_dict['endtime']
    t_start = datetime.datetime(d0.year,d0.month,d0.day,d0.hour,d0.minute,d0.second)
    t_end = datetime.datetime(d1.year,d1.month,d1.day,d1.hour,d1.minute,d1.second)
    # verbatim
    if verbatim > 1:
        print("\t%s read"%file_path)
    return s, fs, t_start, t_end, length_n

def read_fish(file_path, config, verbatim=0):
    """ Read fish data.
    /!\ Signature of this function should not be modified and is similar for all applications (i,e, read_XXX)
    INPUT:
    - file_path: to file containing data (.dat here)
    - config: config dictionnary according to project formating
    - verbatim
    OUTPUT:
    - s: numpy array containing the signal read
    - fs: sampling_rate
    - t_start and t_end: as datetime.datetime objects
    - length_n
    """
    # Read and check reading
    if not isfile(file_path):
        print("No file at",file_path)
        return None, None
    s, fs = soundfile.read(file_path)
    # Get signals metadata
    [_,date,time] = file_path.split('/')[-1].split('__')[0].split('_')
    [y0, m0, d0] = np.array(date.split('-')).astype(int)
    [h0, min0, s0] = np.array(time.split('-')).astype(int)
    length_n,=np.shape(s)
    t_start = datetime.datetime(y0,m0,d0,h0,min0,s0)
    t_end = t_start + datetime.timedelta(seconds=length_n/fs)
    # Correction instrumentale
    SH = -155           # SH en dB ref 1V/µPai, DO NOT CHANGE
    G = 14.7            # gain en dB, DO NOT CHANGE
    D = 2.5             # dynamique de mesure en V, DO NOT CHANGE
    s = s * pow(10, -G/20) * D * pow(10, -SH/20)
    # verbatim
    if verbatim > 1:
        print("\t%s read"%file_path)
    return s, fs, t_start, t_end, length_n

def request_merapi(config, tStart, duration, verbatim=0):
    """ Request Merapi function.
    Gets the signature only, not full signal
    /!\ Signature of this function should not be modified and is similar for all applications (i,e, request_XXX)
    INPUT:
    - config: config dictionnary according to project formating
    - tStart: datetime object
    - duration: in sec
    - verbatim
    OUTPUT:
    - signature: numpy array containing the signal read
    - fs: sampling_rate
    """
    if verbatim > 1:
        debug = True
    else:
        debug = False

    try:
        client = Client(user=config.data_to_analyze['reading_arguments']['user'],
                        host=config.data_to_analyze['reading_arguments']['host'],
                        port=config.data_to_analyze['reading_arguments']['port'],
                        debug=debug)
    except Exception as inst:
        print('Impossible to reach client ')
        print('--', inst)
        return 0, []

    delta_t = eval(config.data_to_analyze['reading_arguments']['delta_t'])
    t = UTCDateTime(tStart)

    try:
        st = client.get_waveforms(  config.data_to_analyze['reading_arguments']['network'],
                                    config.data_to_analyze['reading_arguments']['station'],
                                    config.data_to_analyze['reading_arguments']['location'],
                                    config.data_to_analyze['reading_arguments']['channel'],
                                    t-delta_t, min(t+duration, t+config.data_to_analyze['reading_arguments']['max_duration']))
    except Exception as inst:
        print('Reading not possible for data: ', t, duration )
        print('--', inst)
        return 0, []

    signature = st[0].data
    fs = st[0].stats['sampling_rate']

    if eval(config.data_to_analyze['reading_arguments']['filtering']):
        signature = filter_data(signature, fs, config.data_to_analyze['reading_arguments']['filtering_frequency'])

    return fs, signature

def read_example(file_path, config, verbatim=0):
    """ Read example data.
    /!\ Signature of this function should not be modified and is similar for all applications (i,e, read_XXX)
    INPUT:
    - file_path: to file containing data (.dat here)
    - config: config dictionnary according to project formating
    - verbatim
    OUTPUT:
    - s: numpy array containing the signal read
    - fs: sampling_rate
    - t_start and t_end: as datetime.datetime objects
    - length_n
    """
    # Read and check reading
    if not isfile(file_path):
        print("No file at",file_path)
        return None, None
    s, fs = soundfile.read(file_path)
    s = s[:,0]

    # Get signals metadata
    [date,_,time] = file_path.split('/')[-1].split('.')[0].split('_')
    [y0, m0, d0] = np.array(date.split('-')).astype(int)
    [h0, min0, s0] = np.array(time.split('-')).astype(int)
    length_n,=np.shape(s)
    t_start = datetime.datetime(y0,m0,d0,h0,min0,s0)
    t_end = t_start + datetime.timedelta(seconds=length_n/fs)

    # verbatim
    if verbatim > 1:
        print("\t%s read"%file_path)
        print(np.shape(s), fs, t_start, t_end, length_n)
        print()
    return s, fs, t_start, t_end, length_n



def requestObservation(config, tStart, duration, pathToRecording, verbatim=0):
    """ Request observation function.
    Gets the signature only, not full signal.
    Uses the reading_function to get the full signal.
    - config: config dictionnary according to project formating
    - tStart: datetime object
    - duration: in sec
    - pathToRecording: for use with data reading function. None is analysis_type
    is sparse_realtime
    - verbatim
    OUTPUT:
    - signature: numpy array containing the signal read
    - fs: sampling_rate
    """

    # If reading from online server
    if config.general['analysis_type'] == "sparse_realtime":
        request_function = config.data_to_analyze['reading_function']
        return request_function(config, tStart, duration, verbatim=verbatim)

    # If reading from local recordings
    else:
        # Read full recording
        reading_function = config.data_to_analyze['reading_function']
        print(reading_function, pathToRecording)
        print(config)
        print(verbatim)
        [data, fs, tStartRecording, tEndRecording, length_n] = reading_function(pathToRecording,config,verbatim=verbatim)
        """
        try:
            tStartSignatureInRecording = (tStart - tStartRecording).total_seconds()
            nStartSignatureInRecording = int(tStartSignatureInRecording*fs)
            nEndSignatureInRecording = nStartSignatureInRecording + int(duration*fs)
            signature = data[nStartSignatureInRecording:nEndSignatureInRecording]
        except:
            fs = 1
            signature = np.array([0])
        """
        if verbatim>3:
            print('Fs=',fs,' duration=',duration)
        tStartSignatureInRecording = (tStart - tStartRecording).total_seconds()
        nStartSignatureInRecording = int(tStartSignatureInRecording*fs)
        if not duration or np.isnan(duration): # added by Glenn as length is coming is as NaN, though it wasn't in Paris in July
            duration = np.min([15.0, len(data)/fs])
        nEndSignatureInRecording = nStartSignatureInRecording + int(duration*fs)
        signature = data[nStartSignatureInRecording:nEndSignatureInRecording]        
        return fs, signature
    

def read_montserrat(file_path, config, verbatim, traceID=None):
    """ Read Montserrat data.
    INPUT:
    - file_path: to file containing data (.dat here):
    - verbatim
    OUTPUT:
    - s: numpy array containing the signal read
    - fs: sampling_rate
    - t_start and t_end: as datetime.datetime objects
    - length_n
    """
    
    s = np.array([0])
    fs = 1 
    t_start = None 
    t_end = None 
    length_n = None
    
    # what traceID are we looking for - should figure out how to write this into the config
    if not traceID:
        fptr = open('./AAA-master/MONTSERRAT/current_traceID.txt','r')
        traceID = fptr.read()
        fptr.close()    

    
    # first check for features file
    mseedbase = os.path.basename(file_path)
    features_filepath = os.path.join('features',mseedbase.replace('.mseed', '.%s.features.pkl' % traceID))
   
    
    # Read and check reading
    #file_path = file_path.replace('eventFiles', '../MONTSERRAT/eventFiles')
    if not isfile(file_path):
        return None, None
    else:
        stream = obspy.read(file_path)
    
    if verbatim > 3:
        print('Got %d traces from %s' % (len(stream), file_path))
    
    stream = stream.select(id=traceID)
    L = len(stream)
    if verbatim > 3:
        print('%d traces match %s' % (L, traceID))     
    print(stream)

    if L==1:
        trace = stream.traces[0]
        
        # Get signal and its metadata
        s = trace.data
        s_dict = trace.stats
        # Gets information about recording
        fs = s_dict['sampling_rate'] 
        length_n=s_dict['npts'] # only change from read_ubinas
        if length_n != np.shape(s)[0]:
            print('problem with signals reading dimensions:', length_n, np.shape(s)[0])
        d0 = s_dict['starttime']
        d1 = s_dict['endtime']
        t_start = datetime.datetime(d0.year,d0.month,d0.day,d0.hour,d0.minute,d0.second)
        t_end = datetime.datetime(d1.year,d1.month,d1.day,d1.hour,d1.minute,d1.second)
        
        # file_path is to a pickle file, I think.
        # we want to put that, or the corresponding trace csv file, into a fixed place so that featureFunctions.py can access it.
        tracecsv = file_path.replace('.mseed','.csv')
        outcsv = './AAA-master/MONTSERRAT/precomputed_metrics.csv'
        if isfile(tracecsv):
            tracedf = pd.read_csv(tracecsv)
            # subset for this traceID
            tracedf = tracedf[tracedf['id']==traceID]
            if len(tracedf.index)==1:
                tracedf.to_csv(outcsv)
            elif isfile(outcsv):
                os.system('rm %s' % outcsv)
        else:
            print('Could not find %s' % tracecsv)
            if isfile(outcsv):
                os.system('rm %s' % outcsv)

    else: 
        print('Found %d traces for %s in %s: ' % (L, traceID, file_path) )

    #return s, fs, t_start, t_end, length_n
    return fs, s
