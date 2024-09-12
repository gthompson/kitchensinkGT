import obspy
import pandas as pd
import os
import sys
sys.path.append('../lib')
from SAM import DSAM 

# Load raw seismic data - and set units accordingly
DATA_DIR = os.path.join('..', '..', 'data')
CATALOG_DIR = os.path.join(DATA_DIR, 'catalogs')
halprojdir = '/home/thompsong/Developer/kitchensinkGT/PROJECTS/MVOcatalog/to_dataframes/'
SDS_DIR = os.path.join(DATA_DIR, 'continuous','SDS')
SAM_DIR = os.path.join(DATA_DIR, 'continuous','SAM')
RESPONSE_DIR = os.path.join(DATA_DIR, 'responses')
aefcat = pd.read_pickle(os.path.join(halprojdir, 'MBWH.aef.fixed.pkl'))
badrows = []
dates = []
times = []
for i,row in aefcat.iterrows():
    d=row['date']
    t=row['time']
    try:
        thisdate = obspy.UTCDateTime("%s/%s/%s %s:%s:%s" %(d[0:4], d[5:7], d[8:10], t[0:2], t[3:5], t[6:8])).datetime
        dates.append(thisdate)
    except:
        badrows.append(i)
aefcat.drop(badrows, inplace=True)
aefcat['date']=dates
aefcat.drop(columns=['time'], inplace=True)
print(aefcat)
#aefcat.to_pickle(os.path.join(halprojdir, 'MBWH.aef.fixed.datetime.pkl'))

aefcat.plot(x='date', y='duration')

import numpy as np
aefcat['mag'] = np.log10(aefcat['amp'])
aefcat.loc[aefcat["mag"] < -2.0, "mag"] = -2.0
aefcat.plot.scatter(x='date', y='mag')

mbwhdatcat = pd.read_pickle(os.path.join(halprojdir, 'mbwh_subclass_magnitude_sorted.pkl'))
print(mbwhdatcat)
mbwhdatcat['date'] = pd.to_datetime(mbwhdatcat[['year', 'month', 'day', 'hour', 'minute', 'second']])

print(mbwhdatcat)

#mbwhdatcat.to_pickle(os.path.join(halprojdir,'mbwh_subclass_magnitude_sorted.datetime.pkl'))

aefcat.set_index('date',inplace=True)
mbwhdatcat.set_index('date',inplace=True)
aef_emag_cat = aefcat.merge(mbwhdatcat, left_index=True, right_index=True, how='inner')
print(aef_emag_cat)

aef_emag_cat.reset_index(inplace=True)

ax = aef_emag_cat.plot.scatter(x='date', y='emag')
ax.set_ylim([-2, 5])

# Can i use vsmTools and create a volcanoSeismicCatalog object? Can I plot eventrate by subclass?
import sys
sys.path.append('../lib')
import vsmTools
from obspy.core.event import Event, Catalog, Origin
from obspy.core.event.magnitude import Amplitude
from obspy.core.event.base import CreationInfo, WaveformStreamID
from obspy.core.event.magnitude import Magnitude, StationMagnitude

def mbwhcat2catalog( df ):
    cat = vsmTools.VolcanoSeismicCatalog() 
    cat.classifications = df['subclass_y'].to_list()
    
    cat.starttime = df.iloc[0]['date']
    cat.endtime = df.iloc[-1]['date']
    #cat.magnitudes = df['emag'].to_list()
    #cat.energy = df['eng'].to_list()
    #cat.amplitude = df['amp'].to_list()
    #cat.peakf = df['amp'].to_list()
    this_st=None
    thistrig=None
    for i, row in df.iterrows():
        #print(row)
        
        origin_object = Origin(time=row['date'])
        amplitude_objects = []
        magnitude_objects = []
        stationmag_objects = []
        sta_mags = []

        sta_amp = row['amp']
        if not np.isnan(sta_amp):
            
            amp_obj = Amplitude(generic_amplitude=row['amp']/1e6, \
                                unit='m/s', waveform_id = WaveformStreamID(seed_string='MV.MBWH..SHZ') )
            amplitude_objects.append(amp_obj)
            
        sta_mag = row['emag']
        sta_mags.append(sta_mag)
        stationmag_objects.append(StationMagnitude(mag=sta_mag, mag_type='M'))
        avg_mag = np.nanmean(sta_mags)
        networkmag_object = Magnitude(mag=avg_mag, mag_type='M')
        magnitude_objects.append(networkmag_object)

        info = CreationInfo(author="MVO_energy_magnitude_2000", creation_time=obspy.UTCDateTime())
        this_event = Event(EventType="not reported", creation_info=info, origins=[origin_object], \
                           amplitudes=amplitude_objects, magnitudes=magnitude_objects, station_magnitudes=stationmag_objects)
        cat.addEvent(thistrig, this_st, this_event)
        cat.streams = []
    return cat
    
cat = mbwhcat2catalog( aef_emag_cat )
print(cat)

cat.save(halprojdir, 'mbwh_volcanoSeismicCatalog')
print(os.listdir(halprojdir))

#cat2=vsmTools.load_catalog(halprojdir, 'mbwh_volcanoSeismicCatalog')

import pandas as pd
cat.plot_eventrate(binsize=pd.Timedelta(days=7))



#aef_emag_cat.to_pickle(os.path.join(halprojdir,'mbwh_subclass_aef_emag.pkl'))

print(aef_emag_cat)

df_h = aef_emag_cat[aef_emag_cat['subclass_x']=='h']
ax2=df_h.plot(x='date', y='emag')
ax2.set_ylim([-2, 5])

