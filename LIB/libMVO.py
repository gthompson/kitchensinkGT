#!/usr/bin/env python
import sys
import numpy as np
import os
from glob import glob
from obspy import read_inventory
#LIBpath = os.path.join( os.getenv('HOME'),'src','kitchensinkGT', 'LIB')
#sys.path.append(LIBpath)
from libseisGT import get_seed_band_code
#import metrics


def change_last_sample(tr):
    # For some SEISAN files from Montserrat - possibly from SEISLOG conversion,
    # the last value in the time series was always some absurdly large value
    # So remove the last sample
    tr.data = tr.data[0:-2]

def swap32(i):
    # Change the endianess
    return struct.unpack("<i", struct.pack(">i", i))[0]

def fix_trace_id(st, shortperiod=False):
    # convenience method to wrap correct_nslc
    for tr in st:
        if not 'history' in tr.stats:
            tr.stats['history'] = list() 
        nslc = correct_nslc(tr.id, tr.stats.sampling_rate, shortperiod=shortperiod)
        tr.id = nslc
        if 'deconvolved' in tr.stats.history or 'calibrated' in tr.stats.history:
            if tr.stats.channel[1] == 'D':
                tr.stats['units'] = 'Pa'
            if tr.stats.channel[1] == 'H':
                tr.stats['units'] = 'm/s'    
        else:
            tr.stats['units'] = 'Counts'

def correct_nslc(traceID, Fs, shortperiod=False):
    # Montserrat trace IDs are often bad. return correct trace ID
    oldnet, oldsta, oldloc, oldcha = traceID.split('.')

    net = 'MV'    
    sta = oldsta
    loc = oldloc
    chan = oldcha
    
    if len(chan)>0:
        if chan[0]=='E' or chan[0]=='S':
            shortperiod=True

    if len(chan)>1:
        if chan[1] in 'ZNE':
            chan = chan[0] + 'H' + chan[1]

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
    orientationcode = ''
        
    if chan[0]=='P' or chan[0:2]=='AP':
        instrumentcode ='D' # infrasound/acoustic
        orientationcode = 'F'
        if chan[1].isnumeric(): # e.g. channel like P5
            loc = chan[1]
        
    if instrumentcode == 'H': # seismic
        if len(chan)==2:
            orientationcode = loc
            loc = ''
        elif len(chan)>2:
            orientationcode = chan[2]
    

    chan = bandcode + instrumentcode + orientationcode

    newID = net + "." + sta + "." + loc + "." + chan
    #print(traceID,'->',newID)
    return newID

def inventory_fix_id_mvo(inv):
    inv[0].code='MV'
    net = inv[0].code
    for station in inv[0].stations:
        sta = station.code
        for channel in station.channels:
            chan = channel.code
            if chan[0] in 'ES':
                shortperiod = True
            if chan[0] in 'BH':
                shortperiod = False
            Fs = channel.sample_rate
            nslc = net + '.' + sta + '..' + chan
            nslc = correct_nslc(nslc, Fs, shortperiod=shortperiod)
            net, sta, loc, chan = nslc.split('.')
            channel.code = chan
        station.code = sta
    inv[0].code = net
    return inv

def load_mvo_inventory(tr, CALDIR):
    this_inv = None
    matchcode = None
    if tr.stats.channel[0] in 'ES':
        matchcode = '[ES]'
    elif tr.stats.channel[0] in 'BH':
        matchcode = '[BH]'
    if not matchcode:
        print("Cannot match trace ID %s ",tr.id)
        return this_inv
    xmlfilepattern = os.path.join(CALDIR, "station.MV.%s..%s*%s.xml" % (tr.stats.station, matchcode, tr.stats.channel[2]) )
    xmlfiles = glob(os.path.join(CALDIR, "station.MV.%s..%s*%s.xml" % (tr.stats.station, matchcode, tr.stats.channel[2]) ))
    N = len(xmlfiles)
    if N==1:
        xmlfile = xmlfiles[0]
        print('Correcting %s with %s' % (tr.id, xmlfile))
        this_inv = read_inventory(xmlfile)  
    return this_inv

if __name__ == '__main__':
    pass
