#!/usr/bin/env python
import sys
import numpy as np
sys.path.append('/Users/thompsong/src/kitchensinkGT/LIB')
from libseisGT import get_seed_band_code
#import metrics


def change_last_sample(tr):
# function tr = change_last_sample(tr)
    # For some SEISAN files from Montserrat - possibly from SEISLOG conversion,
    # the last value in the time series was always some absurdly large value
    # This sections was to change last value to be the mean of the rest
    #a = tr.data
    #a = np.delete(a,[np.size(a) - 1])
    #m = np.mean(a)
    #tr.data[-1] = m
    return tr

def swap32(i):
# function tr.data = swap(tr.data)
    return struct.unpack("<i", struct.pack(">i", i))[0]

def fix_trace_id(st, shortperiod=False):
    # convenience method to wrap correct_nslc
    for tr in st:
        nslc = correct_nslc(tr.id, tr.stats.sampling_rate, shortperiod=shortperiod)
        tr.id = nslc

def correct_nslc(traceID, Fs, shortperiod=False):
    # Montserrat trace IDs are often bad. return correct trace ID
    oldnet, oldsta, oldloc, oldcha = traceID.split('.')

    net = 'MV'
    chan = oldcha
    sta = oldsta
    loc = oldloc
    #sta = oldsta.strip()
    #chan = oldcha.strip()
    
    if loc == 'J':
        loc = ''
        
    if loc == '--':
        loc = ''
                
    # channel code is bandcode + instrumentcode + orientationcode
    bandcode = get_seed_band_code(Fs, shortperiod=shortperiod)
    #print(bandcode)
    
    instrumentcode = 'H' # seismic velocity sensor is default
        
    if chan[0]=='P' or chan[0:2]=='AP':
        instrumentcode ='D' # infrasound/acoustic
        orientationcode = 'F'
        if chan[1].isnumeric(): # e.g. channel like P5
            loc = chan[1]
        
    if instrumentcode == 'H': # seismic
        if len(chan)==2:
            orientationcode = loc
            loc = ''
        else:
            orientationcode = chan[2]

    chan = bandcode + instrumentcode + orientationcode

    return net + "." + sta + "." + loc + "." + chan


if __name__ == '__main__':
    pass