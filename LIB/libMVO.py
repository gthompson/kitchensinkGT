#!/usr/bin/env python
import sys
import numpy as np
import os
import pandas as pd
from glob import glob
from obspy import read_inventory
import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as cf
LIBpath = os.path.join( os.getenv('HOME'),'src','kitchensinkGT', 'LIB')
sys.path.append(LIBpath)
from libseisGT import get_seed_band_code


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

def correct_nslc(traceID, Fs, shortperiod=None):
    # Montserrat trace IDs are often bad. return correct trace ID
    
    # special case - based on waveform analysis, this trace is either noise or a copy of MV.MBLG..SHZ
    if traceID == '.MBLG.M.DUM':
        newID = 'MV.MBLG.1.SHZ'
        return newID
    
    
    oldnet, oldsta, oldloc, oldcha = traceID.split('.')

    net = 'MV'    
    sta = oldsta
    loc = oldloc
    chan = oldcha
    
    # channel code is bandcode + instrumentcode + orientationcode
    if len(chan)>0:
        if chan[0]=='E' or chan[0]=='S':
            shortperiod=True
        if chan[0]=='B' or chan[0]=='H':  
            shortperiod=False    
    bandcode = get_seed_band_code(Fs, shortperiod=shortperiod) # not sure here if BB or SP sensor
    instrumentcode = 'H' # 'H' seismic velocity sensor is default
    orientationcode = 'Z'
    if len(chan)>1:
        instrumentcode = chan[1] 
    if len(chan)>2:
        orientationcode = chan[2]    
    
    # Montserrat BB network 1996-2004 had weirdness like
    # BB stations having channels 'SB[Z,N,E]' and
    # SP stations having channels 'S [Z,N,E]'
    # location code was usually 'J' for seismic, 'E' for pressure
    # channel was 'PRS' for pressure
    # there were also 'A N' channels co-located with single-component Integra LA100s, so perhaps those were some other
    # type of seismometer, oriented North?
    # let's handle these directly here
    if len(chan)==2:
        if len(loc)==1:
            chan=chan+loc # now length 3
            orientationcode = loc
            if not loc.isnumeric():
                loc=''   
    if len(chan)==3:
        if chan[0:2]=='SB':
            bandcode = get_seed_band_code(Fs, shortperiod=False) # just because we know it is BB sensor
            instrumentcode = 'H'
        elif chan[0:2]=='S ':
            bandcode = get_seed_band_code(Fs, shortperiod=True) # just because we know it is SP sensor
            instrumentcode = 'H'
        elif chan[0:3]=='PRS':
            bandcode = get_seed_band_code(Fs, shortperiod=True) # just because we know it is SP sensor
            instrumentcode = 'D'
            orientationcode = 'F'
        elif chan[0:3]=='A N':
            bandcode = get_seed_band_code(Fs, shortperiod=True) # just because we know it is SP sensor
            instrumentcode = 'N'
            orientationcode = 'Z'
            #loc = '10'
        elif chan[0]=='P' or chan[0:2]=='AP' or chan=='PRS': # just because we know it is SP sensor
            bandcode = get_seed_band_code(Fs, shortperiod=True)
            instrumentcode ='D' # infrasound/acoustic
            orientationcode = 'F'
            if chan[1].isnumeric(): # e.g. channel like P5
                loc = chan[1]      
        elif chan[1] in 'ZNE': # seismic component in wrong position
            bandcode = get_seed_band_code(Fs, shortperiod=shortperiod)
            instrumentcode = 'H'  
            orientationcode = chan[1]
        elif len(chan)==2 and instrumentcode == 'H': # e.g. .MBLY.5.SH
            orientationcode = loc
            loc = ''  
    if len(loc)>0:
        if loc[0]=='J' or loc[0]=='E' or loc == "--":
            if len(loc)>1:
                loc=loc[1:]
            else:
                loc=''   

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
    if len(tr.stats.channel)<3:
        return this_inv
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

def parse_STATION0HYP(station0hypfile):
    """ file sample
    012345678901234567890123456 
      MWNH1644.53N 6211.45W 407
      MWNL1644.53N 6211.45W 407
      MWZH1644.53N 6211.45W 407
      MWZL1644.53N 6211.45W 407
      0123456789012345678901234 
    """
    station_locations = []
    if os.path.exists(station0hypfile):
        fptr = open(station0hypfile,'r')
        for line in fptr:
            line = line.strip()
            if len(line)==25:
                station = line[0:4]
                latdeg = line[4:6]
                if latdeg != '16':
                    print(line)
                    continue
                latmin = line[6:8]
                latsec = line[9:11]
                hemisphere = line[11]
                if hemisphere == 'N':
                    latsign = 1
                elif hemisphere == 'S':
                    latsign = -1 
                else:
                    #print(line)
                    continue
                londeg = line[13:15]
                lonmin = line[15:17]
                lonsec = line[18:20]
                lonsign = 1
                if line[20]=='W':
                    lonsign = -1
                elev = line[21:25].strip()
                station_dict = {}
                station_dict['name'] = station
                station_dict['lat'] = (float(latdeg) + float(latmin)/60 + float(latsec)/3600) * latsign
                station_dict['lon'] = (float(londeg) + float(lonmin)/60 + float(lonsec)/3600) * lonsign
                station_dict['elev'] = float(elev)
                
                station_locations.append(station_dict)
        fptr.close()
        return pd.DataFrame(station_locations)

#

def add_station_locations(st, station_locationsDF):
    for tr in st:
        df = station_locationsDF[station_locationsDF['name'] == tr.stats.station]
        df = df.reset_index()
        if len(df.index)==1:
            row = df.iloc[0]
            tr.stats['lat'] = row['lat']
            tr.stats['lon'] = row['lon']
            tr.stats['elev'] = row['elev']
        else:
            #print('Found %d matching stations' % len(index))
            continue

#

def plot_station_amplitude_map(st, station0hypfile=None, outfile=None):
    if not station0hypfile:
        return
    
    station_locationsDF = parse_STATION0HYP(station0hypfile)
    add_station_locations(st, station_locationsDF)  
    
    # draw Montserrat coastline
    extent = [-62.27, -62.12, 16.67, 16.82]
    #central_lon = np.mean(extent[:2])
    #central_lat = np.mean(extent[2:])
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1,1,1, projection=crs.PlateCarree())
    ax.set_extent(extent, crs=crs.PlateCarree())
    ax.add_feature(cf.COASTLINE)
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    
    l = len(st)
    amps = np.zeros((l, 1))
    lons = np.zeros((l, 1))
    lats = np.zeros((l, 1))
    #calibs = []
    textlabels = []
    for i, tr in enumerate(st):
        lons[i] = tr.stats.lon
        lats[i] = tr.stats.lat
        amps[i] = tr.stats.metrics.peakamp
        #calibs.append(tr.stats.calib)
        textlabels.append(tr.stats.station)
    sizes = amps/max(amps) 
    #print('calibs=',calibs)
    g = ax.scatter(x=lons, y=lats,color="red",s=sizes*300,alpha=0.8,transform=crs.PlateCarree());
    g.set_facecolor('none');
    g.set_edgecolor('red');    
    for i,thislabel in enumerate(textlabels):
        ax.text(lons[i],lats[i],thislabel,transform=crs.PlateCarree());
            
    # add the volcano
    ax.scatter(-62.1833326, 16.7166638,  300, marker='*', transform=crs.PlateCarree());
         
    # save or show
    if outfile:
        plt.savefig(outfile,  bbox_inches='tight');
    else:
        plt.show();
        


if __name__ == '__main__':
    pass
