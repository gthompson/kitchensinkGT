# Library for converting Well data from CS to pickle files & MiniSEED
import numpy as np
from obspy import UTCDateTime, Stream, Trace
import os
import pandas as pd

# Define phase 2 lookup table & conversions
# From "2022-12-03_ Field Sheet for Deployment of Groundwater Equipment at NASA_Part_II.pdf" 
def get_transducers_dataframe(paths):
    import os
    import pandas as pd
    if os.path.isfile(paths['transducersCSVfile']):
        transducersDF = pd.read_csv(paths['transducersCSVfile'])
    else:
        
        phase2_startdate = UTCDateTime(2022,7,21,14,7,0)
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
        transducersDF.to_csv(paths['transducersCSVfile'])
    return transducersDF

def list_cs_files(TOPDIR, ext='.csv'):
    import glob
    files = []
    uploaddirs = sorted(glob.glob(os.path.join(TOPDIR, '20??????')))
    for uploaddir in uploaddirs:
        #print(uploaddir)
        files_100Hz = sorted(glob.glob(os.path.join(uploaddir, '100hz/*%s' % ext)))
        #print(len(files_100Hz))
        files.extend(files_100Hz)
        files_baro = sorted(glob.glob(os.path.join(uploaddir,  'Baro/*%s' % ext)))
        #print(len(files_baro))
        files.extend(files_baro)
        files_20Hz = sorted(glob.glob(os.path.join(uploaddir,  '20hz/*%s' % ext)))
        #print(len(files_20Hz))
        files.extend(files_20Hz)
    return files    

# Generate complete list of TOB3 files (raw TOB3 files from CS dataloggers)
def list_loggernet_tob3_files(TOB3_DIR):
    return list_cs_files(TOB3_DIR, ext='.dat')

# Generate complete list of LoggerNet CSV files (converted from TOB3 files with LoggerNet)
def list_loggernet_csv_files(TOB3_DIR, ext='.csv'):
    return  list_cs_files(CSVDIR, ext=ext)


def cast_dataframe(dfcsv):
    dfcsv['TIMESTAMP'] = pd.to_datetime(dfcsv.TIMESTAMP)
    dfcsv['RECORD'] = dfcsv['RECORD'].astype(int)
    for col in dfcsv.columns[2:]:
        dfcsv[col] = dfcsv[col].astype(float)
    #return dfcsv

def read_LoggerNet_csv(csvfile):
    dfcsv = pd.read_csv(csvfile, 
                #dtype={'TOA5':str, '100hz_Sensors':int, 'CR6':float, '18084':float, 
                #        'CR6.Std.12.01':float, 'CPU:VWIRE305_100hz.CR6':float, '20853':float, 'DynamicFreq':float}, 
                parse_dates=['TOA5'])
    dfcsv.columns = dfcsv.iloc[0]
    dfcsv=dfcsv.iloc[3:]
    cast_dataframe(dfcsv)
    return dfcsv

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
    return psi * 6894.76

def psi2feet(psi):
    return 2.31 * psi

def psi2inches(psi):
    return 2.31 * psi * 12

def localtime2utc(this_dt):
    hours = 4
    if this_dt>UTCDateTime(2022,11,6,2,0,0):
        hours = 5
    localTimeCorrection = 3600 * hours
    return this_dt + localTimeCorrection
   
'''
def convert2units(st, transducersDF):
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
'''            
def correct_csvfiles(csvfiles, paths, converted_by_LoggerNet=True, MAXFILES=None, keep_existing=True):
    
    if not MAXFILES or MAXFILES>len(csvfiles):
        MAXFILES = len(csvfiles)
    
    transducersDF = get_transducers_dataframe(paths)

    lod = []
    for filenum, rawcsvfile in enumerate(csvfiles[0:MAXFILES]):
        csvbase = os.path.basename(rawcsvfile)
        print('File %d of %d: %s' % ((filenum+1), MAXFILES, rawcsvfile))
        dirname = os.path.basename(os.path.dirname(rawcsvfile))
        uploaddir = os.path.basename(os.path.dirname(os.path.dirname(rawcsvfile)))
        correctedcsvfile = os.path.join("%s.%s.%s" % (os.path.basename(uploaddir), dirname, csvbase))
        if os.path.isfile(correctedcsvfile) & keep_existing:
            print('- Already DONE')
            df2 = pd.read_csv(correctedcsvfile)
            if converted_by_LoggerNet:
                cast_dataframe(df2) # probably not needed for .py_csv files & might even crash
        else:
            print('- Reading')
            try:
                if converted_by_LoggerNet:
                    df2 = read_LoggerNet_csv(rawcsvfile) # definitely not needed for .py_csv files
                else:
                    df2 = pd.read_csv(rawcsvfile)
            except:
                print('Failed to read %s' % rawcsvfile)
                os.rename(rawcsvfile, rawcsvfile + '.bad')
                continue

            print('- Applying calibration equations')
            for col in df2.columns:
                #print(col)
                if isinstance(col,str) and (col[0:2]=='12' or col[0:2]=='21'):
                    this_transducer = transducersDF[(transducersDF['serial']) == col]
                    #print(this_transducer)
                    if len(this_transducer.index)==1:
                        this_transducer = this_transducer.iloc[0].to_dict()
                        #print(this_transducer)
                        df2[col] = compute_psi(df2[col].to_numpy(), this_transducer)
            print('- writing corrected data to %s' % correctedcsvfile)       
            df2.to_csv(correctedcsvfile)

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
        if timediff.seconds>4*60*60: # files should be no more than 4 hours # SCAFFOLD: SHOULD MARK THESE FILES
            print('Problem likely with start time. Filter out all data more than 4 hours before end')
            df2 = df2[df2['TIMESTAMP']>endtime-pd.to_timedelta(4, unit='h')]
            df2 = df2[df2['TIMESTAMP']<=endtime]
            numValidTimes = len(df2.index)
            passed = False

        # check clock drift - this is what all the GPS CSV files are for
        if converted_by_LoggerNet:
            gpscsv = correctedcsvfile.replace('.csv','_gps.csv')
        else:
            gpscsv = correctedcsvfile.replace('.py_csv','_gps.py_csv')
        if not os.path.isfile(gpscsv):
            gpsdf = measureClockDrift(df2)
            if not gpsdf.empty:
                # write out
                gpsdf.to_csv(gpscsv, index=False)
                passed = False
        else:
            gpsdf = pd.read_csv(gpscsv)

        print('- DONE\n')
        
        # Constructing row of lookup table
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
        #print(lookuptableDF)
        #if passed:
        #    print('- writing to SDS')
        #    convert2sds(df2, SDS_TOP)

    lookuptableDF.to_csv(paths['lookuptable'], index=False)   


def removed_unnamed_columns(df):
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    return df

def convert2sds(df, sdsobj, transducersDF, dryrun=False): # I think df here is supposed to be from a single picklefile
    #print('***')
    #print(df.columns)  
    #print('***')  
    local_startt = UTCDateTime(df.iloc[0]['TIMESTAMP'])
    nextt = UTCDateTime(df.iloc[1]['TIMESTAMP'])
    dt = nextt-local_startt
    utc_startt = localtime2utc(local_startt)
    if utc_startt > UTCDateTime():
        return
    print('local ', local_startt, '-> UTC ', utc_startt)
    st = Stream()    
    #print('***')
    #print(df.columns)    
    for col in df.columns[2:]:
        #print('Processing column %s' % col)
        this_transducer = transducersDF[(transducersDF['serial']) == col]
        #print('***')
        #print(this_transducer)
        #print('***')
        if len(this_transducer.index)==1:
            this_transducer = this_transducer.iloc[0].to_dict()
            tr = Trace()
            tr.id = this_transducer['id']
            tr.stats.starttime = utc_startt
            tr.stats.delta = dt  
            tr.data = np.array(df[col])           
            print(f"sampling rate = {tr.stats.sampling_rate}")
            if int(tr.stats.sampling_rate)==20:
                if tr.stats.channel[0]=='H':
                    tr.stats.channel="B%s" % tr.stats.channel[1:]
            if int(tr.stats.sampling_rate)==1:
                if tr.stats.channel[0]=='H' or tr.stats.channel[0]=='B':
                    tr.stats.channel="L%s" % tr.stats.channel[1:]
            #print(tr)
            st.append(tr)
    print('Final Stream object to write')
    print(st)
    if dryrun:
        return False
    else:
        sdsobj.stream = st
        successful = sdsobj.write()
        if successful:
            print("Wrote whole Stream object to SDS")
        else: 
            print("Failed to write entire Stream object:")
            print(df)
        return successful

def remove_overlaps(a):
    c = np.where(np.diff(a) <= 0)[0]
    bad_i = []
    for d in c:
        bool_i = a[0:d] >= a[d+1]
        new_bad_i = list(np.where(bool_i)[0] ) 
        new_bad_i.append(d)
        bad_i = bad_i + new_bad_i
    return list(set(bad_i))


def convert2sds_badfile(df, sdsobj, transducersDF, dryrun=False): # I think df here is supposed to be from a single picklefile
    # NEED TO MODIFY THIS TO DEAL WITH TIME SKIPS
    # No time skip should be more than half a sample.
    # Backward slips are harder to deal with, since time will overlap. Overwrite previous samples when that happens.

    if dryrun:
        mseeddirname = os.path.join(sdsobj.topdir, 'well_mseed')
        if not os.path.isdir(mseeddirname):
            os.makedirs(mseeddirname)

    print('df original length: ', len(df))

    # get rid of anything not within 4 hours of last row
    d = pd.to_datetime(df['TIMESTAMP'])
    df2 = df[ d > d[d.size-1]-pd.Timedelta(hours=4) ]
    #good_i = [i for i in range(len(df2))]
    print('df2 original length: ', len(df2))
    
    # find correct dt using mode
    t = np.array( [UTCDateTime(df2.iloc[i]['TIMESTAMP']).timestamp for i in range(len(df2)) ])
    dt_array = np.diff(t)
    dt = np.median(dt_array)
    #dt_series = pd.Series(dt_array)
    #print(dt_series.describe())
    #dt = dt_series.mode()[0]

    # maybe we only care about times where t slipped backwards?
    # if it goes forward, do we care?
    # do we insert interpolate onto regular sampling?
    # so just generate a list of -ve diffs
    # then process each of those to discover overlaps

    # Now we have expected dt, we need to check ...
    #good_i = range(len(dt_list)+1) # last row always assumed good
    bad_i = remove_overlaps(t)     
    print('bad_i length: ', len(bad_i))
    df2 = df2.drop(df2.index[bad_i]) 
    print('df2 new length: ', len(df2))
    del t

    # now all times should be ascending. Resample.#
    if dt > 0.0099 and dt < 0.0101:
        resampleStr = '10ms'
    elif dt > 0.0499 and dt < 0.0501:
        resampleStr = '50ms'
    elif dt > 0.999 and dt < 1.001:
        resampleStr = '1s'
    else:
        print('delta_t not recognized: ', dt)
        

    df2['datetime'] = pd.to_datetime(df2['TIMESTAMP']) 
    df2.set_index('datetime', inplace=True)

 
    df2_resamp = df2.resample(resampleStr).asfreq()
    df2_resamp.reset_index(drop=True) # datetime index gone

    local_startt = UTCDateTime(df2_resamp.iloc[0]['TIMESTAMP'])
    utc_startt = localtime2utc(local_startt)
    if utc_startt > UTCDateTime():
        return
    print('local ', local_startt, '-> UTC ', utc_startt) 
            
    st = Stream()      
    for col in df2_resamp.columns[2:]:
        #print('Processing column %s' % col)
        this_transducer = transducersDF[(transducersDF['serial']) == col]
        #print('***')
        #print(this_transducer)
        #print('***')
        if len(this_transducer.index)==1:
            this_transducer = this_transducer.iloc[0].to_dict()
            tr = Trace()
            tr.id = this_transducer['id']
            tr.stats.starttime = utc_startt
            tr.stats.delta = dt  
            tr.data = np.array(df[col])           
            #print(f"sampling rate = {tr.stats.sampling_rate}")
            if int(tr.stats.sampling_rate)==20:
                if tr.stats.channel[0]=='H':
                    tr.stats.channel="B%s" % tr.stats.channel[1:]
            if int(tr.stats.sampling_rate)==1:
                if tr.stats.channel[0]=='H' or tr.stats.channel[0]=='B':
                    tr.stats.channel="L%s" % tr.stats.channel[1:]
            #print(tr)
            st.append(tr)
    print('Final Stream object to write')
    print(st)
    
    if dryrun:
        print('dryrun')
        for tr in st:
            mseedfilename = os.path.join(mseeddirname, '%s.%s.%s.ms' % (tr.id, tr.stats.starttime.strftime('%Y%m%d_%H%M%S'), tr.stats.endtime.strftime('%H%M%S') ) )
            print('SCAFFOLD: writing ',mseedfilename)
            tr.write(mseedfilename, format='MSEED') #, encoding=5)
        return True
    else:
        sdsobj.stream = st
        successful = sdsobj.write()
        if successful:
            print("Wrote whole Stream object to SDS")
        else: 
            print("Failed to write entire Stream object:")
            print(df)
        return successful
    

def convert2mseed(df, MSEED_DIR, transducersDF): # I think df here is supposed to be from a single picklefile
    #print('***')
    #print(df.columns)  
    #print('***')
    mseeddirname = os.path.join(MSEED_DIR, 'good')
    if not os.path.isdir(mseeddirname):
        os.makedirs(mseeddirname)  
    local_startt = UTCDateTime(df.iloc[0]['TIMESTAMP'])
    nextt = UTCDateTime(df.iloc[1]['TIMESTAMP'])
    dt = nextt-local_startt
    utc_startt = localtime2utc(local_startt)
    if utc_startt > UTCDateTime():
        return
    print('local ', local_startt, '-> UTC ', utc_startt)
    st = Stream()    
    #print('***')
    #print(df.columns)    
    for col in df.columns[2:]:
        #print('Processing column %s' % col)
        this_transducer = transducersDF[(transducersDF['serial']) == col]
        #print('***')
        #print(this_transducer)
        #print('***')
        if len(this_transducer.index)==1:
            this_transducer = this_transducer.iloc[0].to_dict()
            tr = Trace()
            tr.id = this_transducer['id']
            tr.stats.starttime = utc_startt
            tr.stats.delta = dt  
            tr.data = np.array(df[col])           
            print(f"sampling rate = {tr.stats.sampling_rate}")
            if int(tr.stats.sampling_rate)==20:
                if tr.stats.channel[0]=='H':
                    tr.stats.channel="B%s" % tr.stats.channel[1:]
            if int(tr.stats.sampling_rate)==1:
                if tr.stats.channel[0]=='H' or tr.stats.channel[0]=='B':
                    tr.stats.channel="L%s" % tr.stats.channel[1:]
            #print(tr)
            st.append(tr)
    print('Final Stream object to write')
    print(st)
    
    return stream2miniseed(st, mseeddirname)

def remove_overlaps(a):
    c = np.where(np.diff(a) <= 0)[0]
    bad_i = []
    for d in c:
        bool_i = a[0:d] >= a[d+1]
        new_bad_i = list(np.where(bool_i)[0] ) 
        new_bad_i.append(d)
        bad_i = bad_i + new_bad_i
    return list(set(bad_i))


def convert2mseed_badfile(df, MSEED_DIR, transducersDF): # I think df here is supposed to be from a single picklefile
    # NEED TO MODIFY THIS TO DEAL WITH TIME SKIPS
    # No time skip should be more than half a sample.
    # Backward slips are harder to deal with, since time will overlap. Overwrite previous samples when that happens.

    mseeddirname = os.path.join(MSEED_DIR, 'bad')
    if not os.path.isdir(mseeddirname):
        os.makedirs(mseeddirname)

    print('df original length: ', len(df))

    # get rid of anything not within 4 hours of last row
    d = pd.to_datetime(df['TIMESTAMP'])
    df2 = df[ d > d[d.size-1]-pd.Timedelta(hours=4) ]
    #good_i = [i for i in range(len(df2))]
    print('df2 original length: ', len(df2))
    
    # find correct dt using mode
    t = np.array( [UTCDateTime(df2.iloc[i]['TIMESTAMP']).timestamp for i in range(len(df2)) ])
    dt_array = np.diff(t)
    dt = np.median(dt_array)
    #dt_series = pd.Series(dt_array)
    #print(dt_series.describe())
    #dt = dt_series.mode()[0]

    # maybe we only care about times where t slipped backwards?
    # if it goes forward, do we care?
    # do we insert interpolate onto regular sampling?
    # so just generate a list of -ve diffs
    # then process each of those to discover overlaps

    # Now we have expected dt, we need to check ...
    #good_i = range(len(dt_list)+1) # last row always assumed good
    bad_i = remove_overlaps(t)     
    print('bad_i length: ', len(bad_i))
    df2 = df2.drop(df2.index[bad_i]) 
    print('df2 new length: ', len(df2))
    del t

    # now all times should be ascending. Resample.#
    if dt > 0.0099 and dt < 0.0101:
        resampleStr = '10ms'
    elif dt > 0.0499 and dt < 0.0501:
        resampleStr = '50ms'
    elif dt > 0.999 and dt < 1.001:
        resampleStr = '1s'
    else:
        print('delta_t not recognized: ', dt)
        

    df2['datetime'] = pd.to_datetime(df2['TIMESTAMP']) 
    df2.set_index('datetime', inplace=True)

 
    df2_resamp = df2.resample(resampleStr).asfreq()
    df2_resamp.reset_index(drop=True) # datetime index gone

    local_startt = UTCDateTime(df2_resamp.iloc[0]['TIMESTAMP'])
    utc_startt = localtime2utc(local_startt)
    if utc_startt > UTCDateTime():
        return
    print('local ', local_startt, '-> UTC ', utc_startt) 
            
    st = Stream()      
    for col in df2_resamp.columns[2:]:
        #print('Processing column %s' % col)
        this_transducer = transducersDF[(transducersDF['serial']) == col]
        #print('***')
        #print(this_transducer)
        #print('***')
        if len(this_transducer.index)==1:
            this_transducer = this_transducer.iloc[0].to_dict()
            tr = Trace()
            tr.id = this_transducer['id']
            tr.stats.starttime = utc_startt
            tr.stats.delta = dt  
            tr.data = np.array(df[col])           
            #print(f"sampling rate = {tr.stats.sampling_rate}")
            if int(tr.stats.sampling_rate)==20:
                if tr.stats.channel[0]=='H':
                    tr.stats.channel="B%s" % tr.stats.channel[1:]
            if int(tr.stats.sampling_rate)==1:
                if tr.stats.channel[0]=='H' or tr.stats.channel[0]=='B':
                    tr.stats.channel="L%s" % tr.stats.channel[1:]
            #print(tr)
            st.append(tr)
    print('Final Stream object to write')
    print(st)
    return stream2miniseed(st, mseeddirname)

def stream2miniseed(st, mseeddirname):
    successful = True
    for tr in st:
        mseedfilename = os.path.join(mseeddirname, '%s.%s.%s.ms' % (tr.id, tr.stats.starttime.strftime('%Y%m%d_%H%M%S'), tr.stats.endtime.strftime('%Y%m%d_%H%M%S') ) )
        if os.path.isfile(mseedfilename):
            continue
        print('SCAFFOLD: writing ',mseedfilename)
        try:
            tr.write(mseedfilename, format='MSEED') #, encoding=5)
        except:
            successful = False	
    return successful
 
