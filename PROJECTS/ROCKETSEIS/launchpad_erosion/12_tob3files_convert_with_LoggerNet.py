#!/usr/bin/env python
# coding: utf-8

# # Convert Phase 2 well data
# Locher Environmental recorded data from multiple transducers in two adjacent wells near SLC 39A at Kennedy Space Center (KSC), between March and November, 2022. This was part of a seismo-acoustic erosion pilot experiment. During phase 2 of the experiment, which began on July 21st, 2022, vibrating wire sensors were used and found to give more accurate water levels. These data were captured on Campbell Scientific dataloggers and recorded in TOB3 binary file format. The purpose of this notebook is to:
# - read these files
# - apply calibration equations
# - write to pickle files (Python binary files that serialize the data variables)
# Calibration constants are defined in the transducers dataframe (section 2). These are copied from the file "A Pz Linear Gage Calc_NASA Sensors.xls". Other metadata included in this dataframe come from sections 6 and 7 in the file "2022-12-03_ Field Sheet for Deployment of Groundwater Equipment at NASA_Part_II.pdf".
# 
# Note on *TIME ZONES*: 
# - Local time is used in Campbell Scientific binary files, and files converted to Pickle.
# - UTC is used in MiniSEED files, to match seismo-acoustic data.
# 
# Note on *UNITS*:
# - Pi
# 
# To do:
# - create a lookup table, matching file name, start time, end time, and SEED trace-ID.
# 
# 
# # 1. Imports, path variables, and function definitions

# raw data on April 1, 2022 from 16:10 to 16:40 UTC. Launch from SLC40 at 16:24 UTC, watched from Titusville
import os, sys, glob, obspy
import numpy as np
import pandas as pd

import platform
osname = platform.system()
print(osname)

cwd = os.getcwd()
sys.path.append(os.path.join(cwd, 'campbell'))
#import read_cs_files as campbell

if osname=='Darwin':
    HOME = os.getenv('HOME')
    DROPBOX_TOP = os.path.join(HOME, 'Dropbox')
    SDS_TOP = os.path.join(DROPBOX_TOP, 'DATA', 'SDS')
    WELLDATA_TOP = os.path.join(DROPBOX_TOP, 'DATA', 'KSC', 'KSC_Well_Seismoacoustic_Data/WellData')
    TOB3_DIR = os.path.join(WELLDATA_TOP, 'Uploads')
    PKL_DIR = os.path.join(WELLDATA_TOP, 'Converted')
elif osname=='Linux':
    HOME = os.getenv('HOME')
    DROPBOX_TOP = '/raid/newhome/thompsong/Dropbox'
    SDS_TOP = os.path.join(HOME, 'SDS')
    if not os.path.isdir(SDS_TOP):
        os.mkdir(SDS_TOP)
    WELLDATA_TOP = os.path.join(DROPBOX_TOP, 'DATA', 'KSC', 'KSC_Well_Seismoacoustic_Data/WellData')
    TOB3_DIR = os.path.join(WELLDATA_TOP, 'Uploads')
    PKL_DIR = os.path.join(HOME, 'Converted') 
elif osname=='Windows':
    HOME = 'C:\\Users\\thompsong'
    DROPBOX_TOP = 'D:\\Dropbox'
    SDS_TOP = "D:\\SDS"
    if not os.path.isdir(SDS_TOP):
        os.mkdir(SDS_TOP)
    WELLDATA_TOP = os.path.join(DROPBOX_TOP, 'DATA', 'KSC', 'KSC_Well_Seismoacoustic_Data', 'WellData')
    TOB3_DIR = os.path.join(WELLDATA_TOP, 'Uploads')
    PKL_DIR = 'D:\\Converted' 
lookuptable = os.path.join(PKL_DIR,'lookuptable.csv') 
sys.path.append(os.path.join(HOME, 'Documents', 'GitHub', 'kitchensinkGT', 'LIB'))
from libseisGT import write_stream2sds, read_sds, open_sds, sds_get_nonempty_traceids, sds_percent_available_per_day    
#DROPBOX_PROJECT_DIR = os.path.join(DROPBOX_TOP, 'PROFESSIONAL/RESEARCH/3_Project_Documents/NASAprojects/201602 Rocket Seismology/202010 KSC Launchpad Erosion')
#EVENT_MSEED_DIR = os.path.join(DROPBOX_TOP, 'DATA', 'KSC', 'KSC_Well_Seismoacoustic_Data/SeismoAcousticData/Events')

if not os.path.isdir(PKL_DIR):
    os.mkdir(PKL_DIR)

def measureClockDrift(df):
    passed = True
    starttime=df.iloc[0]['TIMESTAMP']
    endtime=df.iloc[-1]['TIMESTAMP']
    timediff = (endtime - starttime)
    nrows=len(df.index)
    numValidTimes = nrows
    secs = np.array([x.timestamp() for x in df['TIMESTAMP']])
    recs = df['RECORD'].to_numpy().astype('int')
    recsdiff = recs[1:-1]-recs[0:-2]
    #print(recsdiff)
    secsdiff = secs[1:-1]-secs[0:-2]   
    print(secsdiff)
    sample_interval = np.nanmedian(secsdiff)
    gps_records = []
    for i in range(0,len(secsdiff)):
        thisdiff = secsdiff[i]
        if thisdiff >= sample_interval*(recsdiff[i]+0.5) or thisdiff <= sample_interval*(recsdiff[i]-0.5):
            # recsdiff suggests consecutive samples, but thisdiff suggests a strange jump in sample time
            # we assume the GPS clock jumped in and reset. these are the times we save for interpolation
            gps_records.append(df.iloc[i]['RECORD'])
            gps_records.append(df.iloc[i+1]['RECORD'])
            try: # this might not always exist if at end of file
                gps_records.append(df.iloc[i+2]['RECORD']) 
            except:
                pass
    df2 = df[df['RECORD'].isin(gps_records)]    
    print('GPS timed samples')
    print(df2)
    return df2
    #input('Press ENTER to continue')


def compute_psi(dig, d):
    psi = np.zeros((len(dig),1))
    #if np.isnan(d['dig0']):
    #    return psi
    for i in range(len(dig)):
        psi[i] = ((dig[i]-d['dig0']) * d['gf'] + (d['tt']-d['tt0'])*d['tf']+(d['bp0']-d['bp']))
        
    #print(level)
    return psi

def psi2pascals(psi):
    psi2kPa = 6.894757
    Pa = psi * psi2kPa * 1000
    return Pa

def psi2depthmetres(psi):
    psi2kPa = 6.894757
    kPa2mH20 = 0.101974
    mH20 = psi * psi2kPa * kPa2mH20
    return mH20

def localtime2utc(this_dt):
    hours = 4
    if this_dt>obspy.UTCDateTime(2022,11,6,2,0,0):
        hours = 5
    localTimeCorrection = 3600 * hours
    return this_dt + localTimeCorrection
    
def convert2sds(df, SDS_TOP):
    local_startt = obspy.UTCDateTime(df.iloc[0]['TIMESTAMP'])
    nextt = obspy.UTCDateTime(df.iloc[1]['TIMESTAMP'])
    dt = nextt-local_startt
    utc_startt = localtime2utc(local_startt)
    if utc_startt > obspy.UTCDateTime():
        return
        
    st = obspy.Stream()    
        
    for col in df.columns[2:]:
        print('Processing column %s' % col)
        this_transducer = transducersDF[(transducersDF['serial']) == col]
        if len(this_transducer.index)==1:
            this_transducer = this_transducer.iloc[0].to_dict()
            tr = obspy.Trace()
            tr.id = this_transducer['id']
            if tr.stats.sampling_rate==20:
                if tr.stats.channel[0]=='H':
                    tr.stats.channel="B%s" % tr.stats.channel[1:]
            tr.stats.starttime = utc_startt
            tr.stats.delta = dt  
            tr.data = np.array(df[col])
            #print(tr)
            st.append(tr)
    print('Final Stream object to write')
    sdsclient = open_sds(SDS_TOP)
    try:
        write_stream2sds(st, sdsclient)
    except:
        for tr in st:
            new_st= Stream()
            new_st.append(tr)
            try:
                write_stream2sds(new_st, sdsclient)
            except:
                pass
    del sdsclient
    
def convert2units(st):
    for tr in st:
        if tr.stats.network=='FL':
            continue
        try:
            this_transducer = transducersDF[(transducersDF['id']) == tr.id] # not working
        except:
            for i,rows in transducersDF.iterrows():
                if row['id'] == tr.id:
                    this_transducer = row
        if this_transducer['type']=='level':
            tr.data = psi2depthmetres(tr.data)
        elif this_transducer['type']=='pressure':
            tr.data = psi2pascals(tr.data)
       


# # 2. Define phase 2 lookup table & conversions
# 
# (From "2022-12-03_ Field Sheet for Deployment of Groundwater Equipment at NASA_Part_II.pdf")


phase2_startdate = obspy.UTCDateTime(2022,7,21,14,7,0)
transducers = []

# Shallow well (HOF-IW0006S)
transducer1 = {'serial':'AirPressureShallow', 'Fs':100, 'sensor':'barometer','shielding':'none',
               'range_kPa_low':100,'range_kPa_high':100,'media':'air', 'type':'pressure', 
               'model':'Keller 0507.01401.051311.07','set_depth_ft':4.46, 'id':'6S.02374.88.HDH'
              } # serial 237488
transducers.append(transducer1)
transducer2 = {'serial':'1226420', 'Fs':100, 'sensor':'vibrating_wire','shielding':'none',
               'range_kPa_low':70,'range_kPa_high':170,'media':'water', 'type':'level', 
               'model':'Geokon 4500AL','set_depth_ft':3.81,
               'dig0':9751, 'gf':-0.006458, 'tt':21.6, 'tt0':21.3, 'tf':-0.008795, 
               'bp':0.0, 'bp0':14.298, 'id':'6S.12264.20.HDD'
              }
transducers.append(transducer2)
transducer3 = {'serial':'1226423', 'Fs':20, 'sensor':'vibrating_wire','shielding':'foam',
               'range_kPa_low':70,'range_kPa_high':170,'media':'water', 'type':'level', 
               'model':'Geokon 4500AL','set_depth_ft':-5.83,
               'dig0':9605, 'gf':-0.006347, 'tt':21.6, 'tt0':22.2, 'tf':-0.004197, 
               'bp':14.504, 'bp0':14.298, 'id':'6S.12264.23.BDD'
              }
transducers.append(transducer3)
transducer4 = {'serial':'1226419', 'Fs':100, 'sensor':'vibrating_wire','shielding':'foam',
               'range_kPa_low':70,'range_kPa_high':170,'media':'water', 'type':'level', 
               'model':'Geokon 4500AL','set_depth_ft':-6.71,
               'dig0':10040, 'gf':-0.006441, 'tt':21.6, 'tt0':21.1, 'tf':-0.010870, 
               'bp':14.504, 'bp0':14.298, 'id':'6S.12264.19.HDD'
              }
transducers.append(transducer4)
transducer5 = {'serial':'1226421', 'Fs':100, 'sensor':'vibrating_wire','shielding':'none',
               'range_kPa_low':70,'range_kPa_high':170,'media':'water', 'type':'level', 
               'model':'Geokon 4500AL','set_depth_ft':-7.71,
               'dig0':9787, 'gf':-0.006724, 'tt':21.6, 'tt0':21.3, 'tf':-0.001145, 
               'bp':14.504, 'bp0':14.298, 'id':'6S.12264.21.HDD'           
               }
transducers.append(transducer5)

# Intermediate well (HOF-IW00061)
transducer6 = {'serial':'AirPressureDeep', 'Fs':100, 'sensor':'barometer','shielding':'none',
               'range_kPa_low':100,'range_kPa_high':100,'media':'air', 'type':'pressure', 
               'model':'Keller 0507.01401.051311.07','set_depth_ft':4.46, 'id':'6I.0XXXX.XX.HDH'
              }
transducers.append(transducer6)
transducer7 = {'serial':'1226429', 'Fs':100, 'sensor':'vibrating_wire','shielding':'none',
               'range_kPa_low':70,'range_kPa_high':170,'media':'water', 'type':'level', 
               'model':'Geokon 4500AL','set_depth_ft':3.71,
               'dig0':9800, 'gf':-0.006428, 'tt':22.6, 'tt0':21.6, 'tf':-0.002384, 
               'bp':0.0, 'bp0':14.298, 'id':'6I.12264.29.HDD'          
              }
transducers.append(transducer7)
transducer8 = {'serial':'2151692', 'Fs':20, 'sensor':'vibrating_wire','shielding':'foam',
               'range_kPa_low':70,'range_kPa_high':170,'media':'water', 'type':'level', 
               'model':'Geokon 4500AL','set_depth_ft':-9.29,
               'dig0':9459, 'gf':-0.008038, 'tt':22.8, 'tt0':21.8, 'tf':-0.007666, 
               'bp':14.296, 'bp0':14.388, 'id':'6I.21516.92.BDD'
              }
transducers.append(transducer8)
transducer9 = {'serial':'2151691', 'Fs':100, 'sensor':'vibrating_wire','shielding':'foam',
               'range_kPa_low':70,'range_kPa_high':170,'media':'water', 'type':'level', 
               'model':'Geokon 4500AL','set_depth_ft':-18.46,
               'dig0':9414, 'gf':-0.008142, 'tt':22.8, 'tt0':21.5, 'tf':-0.008742, 
               'bp':14.296, 'bp0':14.388, 'id':'6I.21516.91.HDD'
              }
transducers.append(transducer9)
transducer10 = {'serial':'2149882', 'Fs':100, 'sensor':'vibrating_wire','shielding':'none',
               'range_kPa_low':70,'range_kPa_high':170,'media':'water', 'type':'level', 
               'model':'Geokon 4500AL','set_depth_ft':-19.29,
               'dig0':9734, 'gf':-0.008075, 'tt':20.7, 'tt0':21.3, 'tf':-0.000675, 
               'bp':14.602, 'bp0':14.389, 'id':'6I.21498.82.HDD'
               }
transducers.append(transducer10)
transducersDF = pd.DataFrame(transducers)
print(transducersDF)
outfile = os.path.join(WELLDATA_TOP, 'transducer_metadata.csv')
transducersDF.to_csv(outfile)


# convert_program='C:\Program Files (x86)\Campbellsci\LoggerNet\csidft_convert.exe'
# infile = os.path.join(WELLDATA_TOP, 'Uploads','20220826','100hz','100hz_Sensors_100Hz1.dat')
# outfile = infile.replace('dat',"xml")
# outfmt = 'csixml'
# cmd = '\"%s\" \"%s\" \"%s\" %s' % (convert_program, infile, outfile, outfmt)
# print(cmd)
# print(os.path.isfile(infile))
# os.system(cmd)

# # 3. Read raw data files, apply calibration equations, write out to pickle and SDS


# Read CSV file converted using LoggerNet
import pandas as pd
keep_existing = True
allcolumns = []
MAXFILES=99999
lod = []
csvfiles = []

def cast_dataframe(dfcsv):
    dfcsv['TIMESTAMP'] = pd.to_datetime(dfcsv.TIMESTAMP)
    dfcsv['RECORD'] = dfcsv['RECORD'].astype(int)
    for col in dfcsv.columns[2:]:
        dfcsv[col] = dfcsv[col].astype(float)
    #return dfcsv

def read_csv(csvfile):
    dfcsv = pd.read_csv(csvfile, 
                #dtype={'TOA5':str, '100hz_Sensors':int, 'CR6':float, '18084':float, 
                #        'CR6.Std.12.01':float, 'CPU:VWIRE305_100hz.CR6':float, '20853':float, 'DynamicFreq':float}, 
                parse_dates=['TOA5'])
    dfcsv.columns = dfcsv.iloc[0]
    dfcsv=dfcsv.iloc[3:]
    cast_dataframe(dfcsv)
    return dfcsv

uploaddirs = sorted(glob.glob(os.path.join(TOB3_DIR, '20??????')))
for uploaddir in uploaddirs:
    print(uploaddir)
    csvfiles_100Hz = sorted(glob.glob(os.path.join(uploaddir, '100hz/*.csv')))
    print(len(csvfiles_100Hz))
    csvfiles.extend(csvfiles_100Hz)
    csvfiles_baro = sorted(glob.glob(os.path.join(uploaddir,  'Baro/*.csv')))
    print(len(csvfiles_baro))
    csvfiles.extend(csvfiles_baro)
    csvfiles_20Hz = sorted(glob.glob(os.path.join(uploaddir,  '20hz/*.csv')))
    print(len(csvfiles_20Hz))
    csvfiles.extend(csvfiles_20Hz)
#print(csvfiles)
print(len(csvfiles))
maxfiles = min([len(csvfiles), MAXFILES])
for filenum, csvfile in enumerate(csvfiles[0:MAXFILES]):
    csvbase = os.path.basename(csvfile)
    print('File %d of %d: %s' % ((filenum+1), maxfiles, csvfile))
    dirname = os.path.basename(os.path.dirname(csvfile))
    uploaddir = os.path.basename(os.path.dirname(os.path.dirname(csvfile)))
    convertedcsvfile = os.path.join(PKL_DIR, "%s.%s.%s" % (os.path.basename(uploaddir), dirname, csvbase))
    if os.path.isfile(convertedcsvfile) & keep_existing:
        print('- Already DONE')
        df2 = pd.read_csv(convertedcsvfile)
        cast_dataframe(df2)
    else:
        print('- Reading')
        try:
            df2 = read_csv(csvfile)
        except:
            print('Failed to read %s' % csvfile)
            os.rename(csvfile, csvfile+'.bad')
            continue

        print('- Applying calibration equations')
        for col in df2.columns:
            print(col)
            if isinstance(col,str) and (col[0:2]=='12' or col[0:2]=='21'):
                this_transducer = transducersDF[(transducersDF['serial']) == col]
                #print(this_transducer)
                if len(this_transducer.index)==1:
                    this_transducer = this_transducer.iloc[0].to_dict()
                    #print(this_transducer)
                    df2[col] = compute_psi(df2[col].to_numpy(), this_transducer)
        print('- writing calibrated data to %s' % convertedcsvfile)       
        df2.to_csv(convertedcsvfile)

    # check start & end time 
    passed = True
    starttime=df2.iloc[0]['TIMESTAMP']
    endtime=df2.iloc[-1]['TIMESTAMP']
    timediff = (endtime - starttime)
    nrows=len(df2.index)
    numValidTimes = nrows
    secs = np.array([x.timestamp() for x in df2['TIMESTAMP']])
    secsdiff = secs[1:-1]-secs[0:-2]
    sample_interval = np.nanmedian(secsdiff)
    #sample_interval2 = timediff.seconds/(nrows-1)  
    if timediff.seconds>4*60*60: # files should be no more than 4 hours
        print('Problem likely with start time. Filter out all data more than 4 hours before end')
        df2 = df2[df2['TIMESTAMP']>endtime-pd.to_timedelta(4, unit='h')]
        df2 = df2[df2['TIMESTAMP']<=endtime]
        numValidTimes = len(df2.index)
        passed = False

    # check clock drift
    gpscsv = convertedcsvfile.replace('.csv','_gps.csv')
    if not os.path.isfile(gpscsv):
        gpsdf = measureClockDrift(df2)
        if not gpsdf.empty:
            # write out
            gpsdf.to_csv(gpscsv)
            passed = False
    else:
        gpsdf = pd.read_csv(gpscsv)
    

    
    print('- DONE\n')
    thisd = {}
    sourcefileparts = csvfile.split('\\')
    sourcefilerelpath = os.path.join(sourcefileparts[-3],sourcefileparts[-2],sourcefileparts[-1])
    thisd['sourcefile']=sourcefilerelpath
    thisd['outputfile']=convertedcsvfile
    thisd['starttime']=starttime.strftime('%Y/%m/%d %H:%M:%S')
    thisd['endtime']=endtime.strftime('%Y/%m/%d %H:%M:%S')
    thisd['hours']=np.round(timediff.seconds/3600.0,2)
    thisd['npts']=nrows
    thisd['nRECS']=df2.iloc[-1]['RECORD']-df2.iloc[0]['RECORD']+1
    thisd['Fs']=1/sample_interval
    if thisd['nRECS']!=nrows:
        passed=False
    if numValidTimes<nrows:
        passed=False
    thisd['numValidTimes']=numValidTimes
    thisd['numTimeSkips']=len(gpsdf.index)/3
    thisd['passed']=passed
    lod.append(thisd)
    lookuptableDF = pd.DataFrame(lod)
    print(lookuptableDF)
    if passed:
        print('- writing to SDS')
        convert2sds(df2, SDS_TOP)

lookuptableDF.to_csv(lookuptable)    

