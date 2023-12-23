#!/usr/bin/env python
# coding: utf-8

# # Segment SDS Archive based on launch event times and generate event web browser
# 
# Launch times come from 'PilotStudy_KSC_Rocket_Launches.xlsx'
# SDS archives are at SDS_TOP and contain data from wells 6I and 6S and from seismo-acoustic stations
# Segmented event waveform files are saved as MiniSEED to EVENT_WAVEFORMS
# 
# 
from IPython import get_ipython
get_ipython().run_line_magic('run', 'header.ipynb')
HTML_DIR = '/var/www/html/thompsong/KSC_EROSION/EVENTS'
PNG_DIR = os.path.join(HTML_DIR, 'images')
EVENT_WAVEFORMS = os.path.join(paths['outdir'], 'EVENTS') # must exist, and Excel file must be here
csv_launches = os.path.join(paths['outdir'], 'PilotStudy_KSC_Rocket_Launches.csv')
csv_launches_detected = os.path.join(paths['outdir'], 'PilotStudy_KSC_Rocket_Launches_detected.csv')


def event_sds2pkl(launchtime, thisSDSobj, EVENT_WAVEFORMS, pretrig=3600, posttrig=3600, overwrite=False):    
    rawfile = os.path.join(EVENT_WAVEFORMS, '%s_raw.pkl' % launchtime.strftime('%Y%m%dT%H%M%S'))
    if os.path.exists(rawfile) and not overwrite:
        print('%s already exists' % rawfile)
    else:
        print('segmenting %s from SDS' % rawfile) 
        startt = launchtime - pretrig
        endt = launchtime + posttrig
        thisSDSobj.read(startt, endt, speed=1)
        st = thisSDSobj.stream
        st.merge(method=0,fill_value=0)

        if len(st)>0:
            try:
                st.write(rawfile, format='pickle')
            except:
                st2 = Stream()
                for tr in st:
                    try:
                        tr.write('tmp.pkl', 'pickle')
                        st2.append(tr)
                    except:
                        print('Failed:\n',tr)
                print(st2)
                if len(st2)>0:
                    st2.write(rawfile, format='pickle')
                #print('Failed to write raw file')
                #print(st)
                #rawfile = None
        else:
            print('Got no data')
            rawfile = None
    return rawfile

def clean(st):
    for tr in st:
        if tr.stats.network != 'FL':
            continue
        tr.detrend('linear')
        tr.filter('highpass', freq=0.2, corners=2) 
        
def apply_calibration_correction(st):
    # calibration correction

    for tr in st:
        if 'countsPerUnit' in tr.stats:
            continue
        else:
            tr.stats['countsPerUnit'] = 1
            if not 'units' in tr.stats:
                tr.stats['units'] = 'Counts'
            if tr.stats.network[0] =='6': # well data
                if tr.stats.channel[2] == 'D':
                    tr.stats.countsPerUnit = 1/LLE.psi2inches(1) # counts (psi) per inch
                    tr.stats.units = 'inches'
                elif tr.stats.channel[2] == 'H':
                    tr.stats.countsPerUnit = 1/6894.76 # counts (psi) per Pa
                    tr.stats.units = 'Pa'
            elif tr.stats.channel[1]=='D':
                tr.stats.countsPerUnit = 720 # counts/Pa on 1 V FS setting
                if tr.id[:-1] == 'FL.BCHH3.10.HD':
                    if tr.stats.starttime < UTCDateTime(2022,5,26): # Chaparral M25. I had it set to 1 V FS. Should have used 40 V FS. 
                        if tr.id == 'FL.BCHH3.10.HDF':
                            tr.stats.countsPerUnit = 8e5 # counts/Pa
                        else:
                            tr.stats.countsPerUnit = 720 # counts/Pa 
                    else: # Chaparral switched to 40 V FS
                        if tr.id == 'FL.BCHH3.10.HDF':
                            tr.stats.countsPerUnit = 2e4 # counts/Pa
                        else:
                            tr.stats.countsPerUnit = 18 # counts/Pa 
                tr.stats.units = 'Pa'

            elif tr.stats.channel[1]=='H':
                tr.stats.countsPerUnit = 3e2 # counts/(um/s)
                tr.stats.units = 'um/s'
            tr.data = tr.data/tr.stats.countsPerUnit
    
def maxamp(tr):
    return np.max(np.abs(tr.data))

def remove_spikes(st):
    SEISMIC_MAX = 0.1 # m/s
    INFRASOUND_MAX = 3000 # Pa
    FEET_MAX = 21 # feet
    #SEISMIC_MIN = 1e-9
    #INFRASOUND_MIN = 0.01
    
    for tr in st:
        ma = maxamp(tr)
        if tr.stats.units == 'm/s':
            tr.data[tr.data > SEISMIC_MAX] = np.nan
            tr.data[tr.data < -1 * SEISMIC_MAX] = np.nan             
        elif tr.stats.units == 'Pa':
            tr.data[tr.data > INFRASOUND_MAX] = np.nan
            tr.data[tr.data < -1 * INFRASOUND_MAX] = np.nan   
        elif tr.stats.units == 'feet':
            tr.data[tr.data > FEET_MAX] = np.nan
            tr.data[tr.data < -1 * FEET_MAX] = np.nan               

from obspy.signal.trigger import coincidence_trigger
from pprint import pprint
import matplotlib.dates as dates
def detectEvent(st, launchtime):
    trig = coincidence_trigger("recstalta", 3.5, 1, st, 3, sta=2, lta=40)
    best_trig = {}
    best_product = 0
    for this_trig in trig:
        thistime = dates.date2num(this_trig['time'])
        this_product = this_trig['coincidence_sum']*this_trig['duration']
        if this_product > best_product:
            best_trig = this_trig
            best_product = this_product
    pprint(best_trig)
    return best_trig['time']

'''
def add_snr(st, assoctime, threshold=1.5):
    nstime = max([st[0].stats.starttime, assoctime-240])
    netime = min([st[0].stats.endtime, assoctime-60])
    sstime = assoctime
    setime = min([st[0].stats.endtime, assoctime+120])    
    for tr in st:
        tr_noise = tr.copy().trim(starttime=nstime, endtime=netime)
        tr_signal = tr.copy().trim(starttime=sstime, endtime=setime)
        tr.stats['noise'] = np.nanmedian(np.abs(tr_noise.data))
        tr.stats['signal'] = np.nanmedian(np.abs(tr_signal.data))
        tr.stats['snr'] = tr.stats['signal']/tr.stats['noise']
        '''

def group_streams_for_plotting(st):
    groups = {}
    stationsWELL = ['6S', '6I']
    for station in stationsWELL:
        stationStream = st.select(network=station)
        #stationIDS = list(set([tr.id for tr in stationStream]))
        groups[station] = stationStream
    streamSA = st.select(network='FL')
    stationsSA = list(set([tr.stats.station for tr in streamSA]))
    for station in stationsSA:
        stationStream = streamSA.select(station=station)
        #stationIDS = list(set([tr.id for tr in stationStream]))
        groups[station] = stationStream
    #print(groups)
    return groups  




# Read launch data into a DataFrame and generate a list of launch times in Python datetime.datetime format


startover = False # starts with original CSV file again
if os.path.isfile(csv_launches_detected) and startover==False:
    launchesDF = LLE.removed_unnamed_columns(pd.read_csv(csv_launches_detected, index_col=None))
else:
    launchesDF = LLE.removed_unnamed_columns(pd.read_csv(csv_launches, index_col=None))
    dt_tmp = pd.to_datetime(launchesDF['Date'] + ' ' +  launchesDF['Time'])
    launchesDF['Date'] = [pdt.to_pydatetime() for pdt in dt_tmp]
    launchesDF.drop(labels='Time', axis=1, inplace=True)
    del dt_tmp

for thisdir in [EVENT_WAVEFORMS, HTML_DIR, PNG_DIR]:
    if not os.path.isdir(thisdir):
        os.makedirs(thisdir)

if not 'rawfile' in launchesDF.columns:
    launchesDF['rawfile'] = None
if not 'corrected_file' in launchesDF.columns: 
    launchesDF['corrected_file'] = None 
if not 'detection_time' in launchesDF.columns:
    launchesDF['detection_time'] = None 
if not 'short_file' in launchesDF.columns:
    launchesDF['short_file'] = None 
if not 'plotted' in launchesDF.columns:
    launchesDF['plotted'] = False     
launchesDF.to_csv(csv_launches_detected) 


# For each launch, segment raw SDS data to multi-trace MiniSEED file in EVENT_WAVEFORMS directory
thisSDSobj = SDS.SDSobj(paths['SDS_TOP'])
for i, row in launchesDF.iterrows():
    launchTime = UTCDateTime(row['Date'])
    print('Processing launch at %s' % launchTime.strftime('%Y-%m-%d %H:%M:%S'))                        
    if not row['corrected_file']:
        if row['rawfile']:
            rawfile = os.path.join(EVENT_WAVEFORMS, row['rawfile'])
        else:
            rawfile = event_sds2pkl(launchTime, thisSDSobj, EVENT_WAVEFORMS, overwrite=False)
            if rawfile:
                launchesDF.at[i, 'rawfile'] = os.path.basename(rawfile)
            else:
                raise Exception('failed to create %s' % rawfile)
        try:
            st = read(rawfile)    
            #st.merge(method=0, fill_value=0)
        except:
            st = Stream()
        print('%s: %d channels' % (rawfile,len(st)))

        # all these functions safe for well traces too
        clean(st) 
        apply_calibration_correction(st)
        remove_spikes(st)
        
        # write corrected event out
        correctedfile =  os.path.join(EVENT_WAVEFORMS, '%s_long.pkl' % launchTime.strftime('%Y%m%dT%H%M%S'))
        print('Writing %s' % correctedfile)
        try:
            st.write(correctedfile, format='PICKLE') # save 2-hour event waveforms
            launchesDF.at[i, 'corrected_file'] = os.path.basename(correctedfile)
        except:
            pass
del thisSDSobj
launchesDF.to_csv(csv_launches_detected) 


print('Detecting events')            
for i, row in launchesDF.iterrows():
    if row['detection_time'] or not startover:
        continue
    launchTime = UTCDateTime(row['Date'])
    if row['corrected_file']:
        correctedfile =  os.path.join(EVENT_WAVEFORMS, row['corrected_file'])
    else:
        continue
    print('Detecting launch at %s' % launchTime.strftime('%Y-%m-%d %H:%M:%S'))                        

    # subset out the seismo-acoustic traces for detection purposes
    if not os.path.isfile(correctedfile):
        print('File not found: ',correctedfile)
        continue
    st = read(correctedfile)
    SA = st.copy().select(network='FL').trim(starttime=launchTime-100, endtime=launchTime+200)
    if len(SA)==0:
        continue

    assocTime = detectEvent(SA, launchTime)
            
    if abs(assocTime-launchTime)>100:
        assocTime=launchTime
    launchesDF.at[i, 'detection_time'] = assocTime
launchesDF.to_csv(csv_launches_detected) 


# Short files

print('Creating short files')
for i, row in launchesDF.iterrows():                       
    if not row['short_file']:
        launchTime = UTCDateTime(row['Date']) 
        print('- Processing launch at %s' % launchTime.strftime('%Y-%m-%d %H:%M:%S')) 
        if row['corrected_file']:
            correctedfile =  os.path.join(EVENT_WAVEFORMS, row['corrected_file'])
        else:
            continue   
        if not os.path.isfile(correctedfile):
            print('File not found: ',correctedfile)
            continue
        # save 3-minute event waveforms
        st = read(correctedfile)
        if len(st)==0:
            continue

        assocTime = row['detection_time']
        if not assocTime:
            continue
        st_short = st.copy()
        st_short.filter('highpass', freq=0.1, corners=2)
        st_short.trim(starttime=assocTime-30, endtime=assocTime+150)
        print(st_short)
        if len(st_short)>0:
            # write shorter corrected event out
            shortfile =  os.path.join(EVENT_WAVEFORMS, '%s_short.pkl' % launchTime.strftime('%Y%m%dT%H%M%S'))
            print('Writing %s' % shortfile)
            try:
                st_short.write(shortfile, format='PICKLE') # save 2-hour event waveforms
            except:
                pass
            launchesDF.at[i, 'short_file'] = os.path.basename(shortfile)
launchesDF.to_csv(csv_launches_detected)         


# Plots

print('Plotting')
for i, row in launchesDF.iterrows():
    if row['plotted']:
        continue
    launchTime = UTCDateTime(row['Date'])
    print('- Plotting launch at %s' % launchTime.strftime('%Y-%m-%d %H:%M:%S')) 
    for ext in ['long', 'short']:
        if ext=='short':
            pklfile = row['short_file']
        else:    
            pklfile = row['corrected_file']
        if not pklfile:
            continue
        pklfile = os.path.join(EVENT_WAVEFORMS, pklfile)
        if not os.path.isfile(pklfile):
            print('File not found: ',pklfile)
            continue
        st = read(pklfile)
        if len(st)==0:
            continue        
        groups = group_streams_for_plotting(st)
        for station, stream_group in groups.items():
            if len(stream_group)>0:
                pngfile = os.path.join(PNG_DIR, '%s_%s_%s.png' % (launchTime.strftime('%Y%m%dT%H%M%S'), station, ext))
                stream_group.plot(equal_scale=False, outfile=pngfile)

launchesDF.to_csv(csv_launches_detected)          

