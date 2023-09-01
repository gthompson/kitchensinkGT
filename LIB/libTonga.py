#!/usr/bin/env python
# coding: utf-8

# # Background
# 
# Inspired by my analysis of Miami Lakes Rboom data to measure the pressure pulse from the Tonga eruption an beamform back to the source, I began drinking in global barometric and infrasound data from IRIS (and Shake.net). Steve then took my measurements and began collating a table of reduced pressure versus other eruptions. This then manifested an invite from Zhigang Peng to collaborate on a Tonga paper to look at all waves emanating from this eruption.
# 
# This library includes all the functions developed, and related imports
# 
# Glenn Thompson, Jan - Apr 2022



####################### IMPORTS  #####################

import os, sys, obspy
import numpy as np
import math
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import obspy
from obspy import UTCDateTime
from datetime import datetime 
LIBpath = os.path.join( os.getenv('HOME'),'src','kitchensinkGT', 'LIB')
sys.path.append(LIBpath)
from libseisGT import inventory2traceid, get_FDSN_inventory, clean_trace, attach_station_coordinates_from_inventory
from obspy.clients.fdsn import Client
#from obspy.core.util import AttribDict
from obspy.geodetics import locations2degrees, degrees2kilometers, kilometers2degrees
from scipy.stats import linregress

####################### CONSTANTS #####################
DATA_ROOT = os.path.join(os.getenv('HOME'), 'DATA', 'HungaTonga')
if not os.path.isdir(DATA_ROOT):
    os.makedirs(DATA_ROOT)
circum_earth_km = 40075
circum_earth_m = circum_earth_km*1000
m_per_deg = circum_earth_m/360.0

# location of Hunga-Tonga Hunga-Ha'apai - used for lightning
sourcelat = -20.536
sourcelon = -175.382

# earthquake origin from USGS
olat = -(20 + 34/60 + 12/3600) 
olon = -(175 + 22/60 + 48/3600)
otime = obspy.core.UTCDateTime('2022-01-15T04:14:45.000000Z') # main eruption time - on day 2


####################### FUNCTIONS FOR SEISMO-ACOUSTICS #####################
def get_fdsn_identifier(fdsnURL):
    prepend=''
    if 'iris' in fdsnURL:
        prepend='iris_'
    if 'shake' in fdsnURL:
        prepend='rboom_'
    if 'geonet' in fdsnURL:
        prepend='nz_'  
    return prepend


# Use get_inventory and get_stream in FDSNtools here
""" def get_inventory_and_waveforms(fdsnURL, searchRadiusDeg, olat, olon, startt, endt, chanstring, data_root, overwrite=False, network='*', station='*', location='*'):
    fdsnClient = Client(fdsnURL)  
    if not os.path.isdir(data_root):
        print('No such directory %s' % data_root)
        return None, None
    
    # STEP 1 - FIND STATIONS BY REQUESTING AN INVENTORY
    prepend = get_fdsn_identifier(fdsnURL)
    stationXmlFile = os.path.join(data_root, 'STATIONXML', '%s%s_within_%.0f_degrees.xml' % (prepend, chanstring, searchRadiusDeg))
    mseedfile = stationXmlFile.replace('STATIONXML', 'MINISEED').replace('.xml','.mseed')
    print(stationXmlFile, mseedfile)
    
    if os.path.exists(stationXmlFile) and overwrite==False:
        inv = obspy.read_inventory(stationXmlFile)
    else:
        inv = get_FDSN_inventory(fdsnClient, startt, stationXmlFile, network, olat, olon,
             searchRadiusDeg, 0, endt - startt, station = station, channel = chanstring ) # get all low rate, outside barometers
    trace_ids = inventory2traceid(inv)
    #print(trace_ids)
    
    # STEP 2 - LOAD CORRESPONDING WAVEFORM DATA
    if os.path.exists(mseedfile) and overwrite==False:
        st = obspy.core.read(mseedfile)
    else:
        st = obspy.core.Stream()
        for trace_id in trace_ids:
            try:
                net, sta, loc, cha = trace_id.split('.')
            except:
                net, sta, cha = trace_id.split('.')
                loc = '*'
            try:
                st0 = fdsnClient.get_waveforms(net, sta, loc, cha, startt, endt)
                print('waveform downloaded for ', trace_id)
                if len(st0)==1:
                    st.append(st0[0])
                else:
                    print('More than 1 Trace')
                    last_id = ''
                    for tr0 in st0:
                        if tr0.id != last_id:
                            st.append(tr0)
                        else:
                            st.remove(st[-1])
            except:
                print('failed to download waveform for ', trace_id)
            
        # save waveform data
        if len(st)>0:
            st.write(mseedfile)
        
    return st, inv """



def smart_merge(st): # a more thorough version exists in libseisGT
    st2 = obspy.core.Stream()
    for tr_original in st:
        tr = tr_original.copy()
        st2 = st2.append(tr)
        try:
            st2.merge(method=1, fill_value=0)
        except:
            st2.remove(tr)
            pass
    return st2



def reconstitute_stream(st, inv, fmin=0.0001):
    T=1.0/fmin
    
    # deal with multiple traces with same trace ID, by appending and merging one at a time,
    # and removing any time there is a fail.
    st2 = smart_merge(st)
        
    # reconstitute as many traces as possible, and return them
    reconstituted = obspy.core.Stream()
    failed_ids = []
    for tr_original in st2:
        tr = tr_original.copy()
        
        # remove gappy Trace objects, filled with 0 above (threshold 10%)
        n = np.count_nonzero(tr.data==0)
        if n/len(tr.data) > 0.1:
            st2.remove(tr_original)
            continue
            
        # high pass filter    
        tr.detrend('linear')
        duration = tr.stats.endtime - tr.stats.starttime
        tr.taper(T/duration)
        tr.stats['maxamp_raw'] = np.max(np.abs(tr.data))
        #tr.filter("highpass", freq=fmin, corners=2)
        pre_filt = [fmin, fmin*2, tr.stats.sampling_rate*0.3, tr.stats.sampling_rate*0.4]
        
        # remove response
        try:
            if tr.stats.channel[1]=='H':
                tr.remove_response(inventory=inv, pre_filt=pre_filt, output='VEL') # was 'DISP'
                tr.stats['units']='m/s'
            else:
                tr.remove_response(inventory=inv, pre_filt=pre_filt )
                tr.stats['units']='Pa'
            reconstituted.append(tr)
            #print(tr.id, ' reconstituted')
        except:
            #print(tr.id, ' NOT reconstituted')
            failed_ids.append(tr.id)
            pass
        tr.detrend('linear')
        tr.stats['maxamp_corrected'] = np.max(np.abs(tr.data))
    if len(failed_ids)>0:
        print('Failed to reconstitute: ', failed_ids)

    #reconstituted.plot(equal_scale=True);
    
    return reconstituted


# FOLLOWING FUNCTIONS NOW IN INVENTORYTOOLS
""" def attach_station_coordinates_from_inventory(inventory, st):
    """ attach_station_coordinates_from_inventory """
    for tr in st:
        for netw in inventory.networks:
            for sta in netw.stations:
                if tr.stats.station == sta.code and netw.code == tr.stats.network:
                    for cha in sta.channels:
                        if tr.stats.location == cha.location_code:
                            tr.stats.coordinates = obspy.core.util.AttribDict({
                                'latitude':cha.latitude,
                                'longitude':cha.longitude,
                                'elevation':cha.elevation})
    return


                            
def attach_distance_to_stream(st, olat, olon):
    for tr in st:
        try:
            alat = tr.stats['coordinates']['latitude']
            alon = tr.stats['coordinates']['longitude']
            distdeg = locations2degrees(olat, olon, alat, alon)
            distkm = degrees2kilometers(distdeg)
            tr.stats['distance'] =  distkm * 1000
        except:
            print('cannot compute distance for %s' % tr.id)
    return """
       

     
def remove_outliers(a, t, f=30, amax=2000, amin=1):
    bad_indices = []
    
    # subset to range amax to amin
    for i, v in enumerate(a):
        if v > amax:
            bad_indices.append(i)
        if v < amin:
            bad_indices.append(i)            
    #print(bad_indices)
    if len(bad_indices)>0:
        for i in bad_indices:
            a[i]=np.NaN 
            t[i]=np.NaN
            
    # subset further based on median
    m = np.median(a)
    bad_indices = []
    for i, v in enumerate(a):
        if v > m * f:
            bad_indices.append(i)
        if v * f < m:
            bad_indices.append(i)            
    #print(bad_indices)
    if len(bad_indices)>0:
        for i in bad_indices:
            a[i]=np.NaN 
            t[i]=np.NaN
    #return a, t
    return
  

  
def plot_amplitude_versus_distance(st, units, reftime=None, cmin=300, cmax=420, duration=7200):
    r = []
    a = []
    t = []
    for tr in st:
        this_r = tr.stats.distance
        if reftime:
            
            Pmax = []
            tmax = []
            stime = reftime + this_r/cmax
            etime = reftime + this_r/cmin + duration
            #print(stime, etime)
            if etime < tr.stats.endtime:
                tr0=tr.copy().trim(starttime=stime, endtime=etime)
                this_a = np.max(np.abs(tr0.data))
                this_t = np.argmax(np.abs(tr0.data))*tr0.stats.delta
                Pmax.append(this_a)
                tmax.append(this_t)
                r.append(this_r/1000)
                a.append(this_a)
                tr.stats['Pmax']=Pmax
                tr.stats['tmax']=tmax
                t.append(this_t) 
                
            this_r = circum_earth_m - this_r
            stime = reftime + this_r/cmax
            etime = reftime + this_r/cmin + duration
            if etime < tr.stats.endtime:
                tr0=tr.copy().trim(starttime=stime, endtime=etime)
                this_a = np.max(np.abs(tr0.data))
                this_t = np.argmax(np.abs(tr0.data))*tr0.stats.delta
                Pmax.append(this_a)
                tmax.append(this_t)
                a.append(this_a)
                tr.stats['Pmax']=Pmax
                tr.stats['tmax']=tmax
                t.append(this_t)    
                r.append(this_r/1000)
        else:
            Pmax = np.max(np.abs(tr.data))
            tmax = np.argmax(np.abs(tr.data))*tr.stats.delta
            a.append(Pmax)
            tr.stats['Pmax']=Pmax
            tr.stats['tmax']=tmax
            t.append(tmax)
            r.append(this_r/1000)
            
    #remove_outliers(a, t) 
    r2 = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
    r2.extend(r)
    PR = 12000
    y = PR/np.power(r2,0.5)
    
    
    # wavefront circumference does not increase like 1/r as assumed because the earth is curved, not flat.
    # so there is a trigonemtric correction to apply. figure this out and use in a better paper, and then
    # estimate attenuation from full energy of N-wave. remove this plot for nopw.
    y2 = []
    r3 = []
    radius_earth = circum_earth_m/(2*math.pi)
    for this_r in r2: 
        real_r = radius_earth * np.sin(this_r/radius_earth)
        #print(this_r, real_r)
        r3.append(real_r)
        y2.append(PR / np.power(real_r, 0.5))
    y2 = np.array(y2)
    r2 = np.array(r2)
    r3 = np.array(r3)
    #print(y2)
    
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(1,1,1)
    ax.loglog(r, a, 'o')
    #ax.loglog(r2, y, 'r:') # this is the plot for 1/sqrt(r)
    #ax.loglog(r2, y2, 'g:')
    ax.set_ylim(10,1000000)
    ax.set_xlim(1,circum_earth_km/4)
    ax.set_ylabel('Maximum amplitude (%s)' % units)
    ax.set_xlabel('Distance (km)')
    plt.grid(which='both')
    #plt.savefig('reduced_pressure.eps')
    
    

    # plot pressure vs distance from model
    r_msvf = 757.0 # km
    P_msvf = 780.0 # Pa
    #Preduced = P_msvf / (1/r_msvf + 1/np.sqrt(crossover_km * r_msvf)) 
 
    r_list, P_list, v_list, t_list, M_list = pressure_vs_distance(r_msvf, P_msvf, crossover_km=100, vambient=340)
    ax.loglog(r_list, P_list, 'k-', linewidth=3)
    print(r_list[:3],P_list[:3])
    r_list, P_list, v_list, t_list, M_list = pressure_vs_distance(r_msvf/3, P_msvf, crossover_km=1, vambient=340)
    ax.loglog(r_list, P_list, 'k--', linewidth=3)    
    print(r_list[:3],P_list[:3])
    # straight line fit
    #fig = plt.figure(figsize=(12,8))
    #ax = fig.add_subplot(1,1,1)
    
  
    r_list = [r for (P,r) in zip(P_list,r_list) if math.isnan(P)==False]
    P_list = [P for P in P_list if math.isnan(P)==False]
    
    r_list = r3[10:]
    P_list = a
    #print(len(r_list), len(P_list))
    
    #print(a)
            
  
    logP = np.log10(np.array(P_list))
    logR = np.log10(np.array(r_list))
    
    logR = logR[:5]
    logP = logP[:5]
    print(logR, logP)
    
    
    slope, logp_intercept, r_value, p_value, slope_std_err = linregress(logR, logP)
    print(slope, logp_intercept)

    logr_intercept = -logp_intercept/slope
    logr2 = np.arange(0.0, np.log10(circum_earth_km), np.log10(circum_earth_km)/1000 )# 1000 steps

    #print(logr2)
    logp2 = slope * logr2 + logp_intercept
    #print(logp2)
    #ax.plot(logR, logP, 'bx') # line
    r2 = np.power(10, logr2)
    p2 = np.power(10, logp2)
    r4 = radius_earth * np.arcsin(r2/radius_earth)
    #ax.loglog(r4, p2, '-')
    ax.set_xlim(1,circum_earth_km/4)
    return
  

  
def plot_reduced_pressure(st, units):
    r = []
    a = []
    for tr in st:
        if 'Pmax' in tr.stats:
            r.append(tr.stats.distance/1000)
            PRmax=tr.stats.Pmax* np.sqrt(tr.stats.distance/1000)
            a.append(PRmax)
            tr.stats['PRmax']=PRmax
    #remove_outliers(a, t, amax=50000, amin=1000)    
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(1,1,1)
    ax.scatter(r, a)
    ax.set_ylabel('Reduced Pressure (%s)' % units)
    ax.set_xlabel('Distance (km)')
    plt.grid()
    plt.show()
    return
    


def yes_or_no(question, default_answer='y', auto=False):
    #while "the answer is invalid":
    if auto==True:
        reply = ''
    else:
        reply = str(input(question+' (y/n): ')).lower().strip()
    if len(reply) == 0:
        reply = default_answer
    if reply[0] == 'y':
        return True
    if reply[0] == 'n':
        return False
   

     
def manually_select_good_traces(st2, trim_mins=0, record_section=False, auto=False, wlen=7200, rmax=999999999):
    r = get_distance_vector(st2)
    st = order_traces_by_distance(st2, r) # sort traces by distance
    st = smart_merge(st)
    st.detrend('linear')
    last_distance = -9999999
    last_max_amp = 999999999999
    stkeep = obspy.core.Stream()
    for tr in st:
        startt = tr.stats.starttime + trim_mins * 60
        endt = tr.stats.endtime - trim_mins + 60 
        if 'arrival' in tr.stats:
            print(tr.stats.arrival)
            startt = tr.stats.arrival['arrivaltime'] - wlen/10
            endt = tr.stats.arrival['arrivaltime'] + wlen + wlen/10
        if endt > startt:
            tr.trim(starttime=startt, endtime=endt)
            tr.plot()
            
            # only recommend to keep a trace if amplitude less than thrice previous
            this_max_amp = np.max(np.abs(tr.data))
            print('max amp = ', this_max_amp)
            this_distance = tr.stats.distance
            if this_distance > rmax: # skip if greater than rmax
                continue
            print('distance = ', this_distance)
            recommend = 'n'
            if this_max_amp < 3 * last_max_amp:
                if record_section: # if record section, only recommend to keep if at least 100 km from previous kept trace
                    if this_distance > last_distance + 200000:
                        recommend = 'y'
                else:
                    recommend = 'y'

            if yes_or_no('keep this trace? (%s)' % recommend, default_answer=recommend, auto=auto):
                stkeep.append(tr.copy())
                last_distance = this_distance
                last_max_amp = this_max_amp
    return stkeep



def select_previously_manually_selected_traces(st, good_ids):
    stkeep = obspy.core.Stream()
    for tr in st:
        if tr.id in good_ids:
            stkeep.append(tr.copy())
    return stkeep



def analyze_clientchan(fdsnClient, chanstring, startt, endt, fmin=0.001, network='*', station='*', location='*', searchRadiusDeg=180.0):
    T = np.round(1/fmin,0)
    
    # raw stream and inventory - note that stationXML and Miniseed files are written if they do not exist already
    this_stream, this_inventory = get_inventory_and_waveforms(fdsnClient, searchRadiusDeg, olat, olon, startt, endt, chanstring, DATA_ROOT, network=network, station=station, location=location, overwrite=False)
    attach_station_coordinates_from_inventory(this_inventory, this_stream)
    attach_distance_to_stream(this_stream, olat, olon)
    
    # reconstituted Stream
    prepend = get_fdsn_identifier(fdsnClient)
    reconfile = os.path.join(DATA_ROOT, 'MINISEED', '%s%s_reconstituted_%.0f.mseed' % (prepend, chanstring, T))  
    this_reconstituted = obspy.Stream()
    if os.path.exists(reconfile):
        choice = input('%s exists. Load this (y/n)? Otherwise data will be reconstituted and saved over this' % reconfile)
        if choice[0].lower() == 'y':
            this_reconstituted = obspy.read(reconfile)
    if len(this_reconstituted)==0:
        this_reconstituted = reconstitute_stream(this_stream, this_inventory, fmin=fmin)
        this_reconstituted.write(reconfile)

    # good traces
    good_ids_pickle = os.path.join(DATA_ROOT, 'PICKLE', '%sgood_%s_ids_within_%.0f_degrees.pkl' % (prepend, chanstring, searchRadiusDeg))
    goodmseedfile = os.path.join(DATA_ROOT, 'MINISEED', '%s%s_good.mseed' % (prepend, chanstring))   
    if os.path.exists(good_ids_pickle):
        with open(good_ids_pickle, 'rb') as f:
            good_ids = pickle.load(f)        
        this_good = select_previously_manually_selected_traces(this_reconstituted, good_ids)
    else:
        this_good = obspy.Stream()
        if os.path.exists(goodmseedfile):
            choice = input('%s exists. Load this (y/n)? Otherwise you have to manually pick traces and overwrite existing file' % goodmseedfile)
            if choice[0].lower() == 'y':
                this_good = obspy.read(goodmseedfile)
        if len(this_good)==0:
            this_good = manually_select_good_traces(this_reconstituted)
            good_ids = []
            for tr in this_good:
                good_ids.append(tr.id)
            with open(good_ids_pickle, 'wb') as f:
                pickle.dump(good_ids, f)    
            this_good.write(goodmseedfile)   
    print('Good traces saved as %s ' % goodmseedfile)
 
    if chanstring[1]=='D':
        units = 'Pa'
        #plot_amplitude_versus_distance(this_good, units)
        #plot_reduced_pressure(this_good, units)
    elif chanstring[1]=='H':
        units = 'm/s'
    #    plot_amplitude_versus_distance(this_good, units)       
    
    good_inv = subset_inv(this_inventory, this_stream, this_good)
    plot_inv(good_inv)
    
    return this_good, good_inv, this_stream, this_inventory



def subset_inv(inv, st, st_subset):
    # subset an inventory based on a stream object which is a subset of another
    try:
        inv_new = inv.copy() # got an error here once that Inventory has no copy(), but it does
        for tr in st:
            if len(st_subset.select(id=tr.id))==0:
                inv_new = inv_new.remove(network=tr.stats.network, station=tr.stats.station, location=tr.stats.location, channel=tr.stats.channel)
        return inv_new
    except:
        print('Failed to subset inventory. Returning unchanged')
        return inv



def plot_inv(inv):
    inv.plot(water_fill_color=[0.0, 0.5, 0.8], continent_fill_color=[0.1, 0.6, 0.1], size=30);
    return
    


def medfilt(x, k):
    """Apply a length-k median filter to a 1D array x.
    Boundaries are extended by repeating endpoints.
    """
    assert k % 2 == 1, "Median filter length must be odd."
    assert x.ndim == 1, "Input must be one-dimensional."
    k2 = (k - 1) // 2
    y = np.zeros ((len (x), k), dtype=x.dtype)
    y[:,k2] = x
    for i in range (k2):
        j = k2 - i
        y[j:,i] = x[:-j]
        y[:j,i] = x[0]
        y[:-j,-(i+1)] = x[j:]
        y[-j:,-(i+1)] = x[-1]
    return np.median(y, axis=1)



def median_downsample_to_1Hz(st, winlec_sec=1.0):
    for tr in st:
        y = medfilt(tr.data, np.round(winlen_secs / tr.stats.delta,0) )
        fs = np.round(tr.stats.sampling_rate,0)
        if fs>=2.0:
            fs_new = tr.stats.sampling_rate/fs
            tr.data = y[0::fs]
            tr.stats.sampling_rate=fs_new
    return

       
     
def plots(this_good, inv, olat, olon, otime, startt=0, endt=0, minf=0, maxf=0, auto=False, cmin=300, cmax=420):  

    st = this_good.copy()
    if endt>startt:
        st.trim(starttime=startt, endtime=endt)
    if minf>0:
        st.filter('highpass', freq=minf, corners=2 )
    if maxf>0:
        st.filter('lowpass', freq=maxf, corners=2 )
    plot_amplitude_versus_distance(st, 'Pa', reftime=otime, cmin=300, cmax=420)
    plot_reduced_pressure(st, 'Pa')   
    plot_inv(inv)

    ntraces = len(st)
    nsamples = len(st[0].data)    
    if nsamples * ntraces > 3600 * 30 * 50: # no more than 1 sample per second for 30 hours for 50 traces
        st.decimate(8)
    print('Select good traces for record section')
    st2 = manually_select_good_traces(st, record_section=True, auto=auto)    
    st2.plot(type='section', scale=10, norm_method='stream', orientation='horizontal', dist_degree=True, ev_coord=[olat,olon]) #, reftime=otime)
    Return



def get_distance_vector(st):
    r = []
    for tr in st:
        r.append(tr.stats.distance)
    return r
     
   
    
def order_traces_by_distance(st, r=[], assert_channel_order=False): 
    st2 = obspy.core.Stream()
    if not r:
        r = get_distance_vector(st)
    if assert_channel_order: # modifies r to order channels by (HH)ZNE and then HD(F123) etc.
        for i, tr in enumerate(st):
            c1 = int(tr.stats.location)/1000000
            numbers = 'ZNEF0123456789'
            c2 = numbers.find(tr.stats.channel[2])/1000000000
            r[i] += c1 + c2
    indices = np.argsort(r)
    for i in indices:
        tr = st[i].copy()
        st2.append(tr)

    return st2

def order_traces_by_id(st):
    ids = []
    for tr in st:
        ids.append(tr.id)
    ids.sort()
    st2 = obspy.core.Stream()
    for id in ids:
        for tr in st:
            if tr.id == id:
                st2.append(tr.copy())
    return st2            


def pick_pressure_onset(st, otime, winsecsmax=7200, speed=310, manual=True):
    print('Pick onset times with left mouse button. Middle button to skip trace.')
    
    for tr in st:
        good_trace = False
        if manual == True:
            winsecs = winsecsmax
            midtime = otime + tr.stats.distance/speed
            mintime = midtime - winsecs/2
            maxtime = midtime + winsecs/2
        
            
        
            while winsecs > 30:
                tr2 = tr.copy().trim(starttime=mintime, endtime=maxtime)
                t = np.arange(0,tr2.stats.delta*len(tr2.data),tr2.stats.delta)
                plt.plot(t, tr2.data)
                plt.title(tr.id)
                plt.show()
                pos = plt.ginput(1, timeout=12.0, show_clicks=True )
                plt.close()
                if len(pos)>0:
                    if len(pos[0])>0:
                        good_trace = True
                        x = pos[0][0]
                        atime = tr2.stats.starttime + x
                        #print(atime)
                        winsecs = winsecs / 5
                        mintime = atime - winsecs/2
                        maxtime = atime + winsecs/2 
                    else:
                        break
                else:
                    break
        
        else:
            tr0=tr.copy().trim(starttime=stime, endtime=etime)
            #Pmax = np.max(np.abs(tr0.data))
            traveltime = np.argmax(np.abs(tr0.data))*tr0.stats.delta
            atime = stime + traveltime
            good_trace = True  
            
            
        if good_trace:
            traveltime = atime - otime
            tr.stats['arrival'] = obspy.core.util.AttribDict({
                    'arrivaltime':atime,
                    'c':tr.stats.distance/traveltime,
                    'traveltime':traveltime})
                #print(tr.stats.arrival)
    return  
        


def regress_arrival_times(st, reftime, outfile=None):
    fh = plt.figure(figsize=(16,16))
    ah = fh.add_subplot(1,1,1)
    atimes=[]
    r=[]
    for i, tr in enumerate(st):   
        if 'arrival' in tr.stats:
            r.append(tr.stats.distance) 
            atimes.append(tr.stats.arrival.arrivaltime - reftime)  
    atimes = np.array(atimes)
    r = np.array(r)
    #m, c = np.polyfit(atimes, r, deg=1)
    wavespeed, d_intercept, r_value, p_value, std_err = linregress(atimes, r)
    ah.plot(atimes, r, 'bx') # data
    t_intercept = -d_intercept/wavespeed
    t = np.arange(min([t_intercept,0.0]), np.max(atimes), np.max(atimes)/1000) # 1000 steps
    d = wavespeed * t + d_intercept
    ah.plot(t, d) # line
    new_otime = reftime + t_intercept
    #print(wavespeed, d_intercept, t_intercept, r_value, p_value, std_err)
    print('infrasonic wave speed = %.1f m/s, origin time = %s' % (wavespeed, new_otime.strftime('%Y-%m-%d %H:%M:%S')))
    ah.set_ylabel('Distance (m)')
    ah.set_xlabel('Time (seconds) after %s' % reftime.strftime('%Y-%m-%d %H:%M:%S'))
    #ah.set_ylim(-1e5, 1e5)
    #ah.set_xlim(-1e3, 3e1)
    ah.grid()
    speed_to_nearest_station = st[0].stats.distance / (st[0].stats.arrival.arrivaltime - new_otime)
    print('Average speed to station %s at %.1f km is %.1f m/s' % (st[0].stats.station, st[0].stats.distance/1000, speed_to_nearest_station ))
    print('To fit earthquake time, speed is %.1f m/s' % (st[0].stats.distance / (st[0].stats.arrival.arrivaltime - reftime)))

    if outfile:
        fh.savefig(outfile)
    return wavespeed, new_otime, t, d



def plot_record_section(st2, reftime, outfile=None, starttime=0, endtime=0, km_min=0, km_max=40075.0/2, slopes = [], scale_factor=1.0, plot_arrivals=False, figsize=(16,16), min_spacing_km=0, reduced_speed=None, normalize=False ):
    r = np.array(get_distance_vector(st2))
    st = order_traces_by_distance(st2, r)


        
    if reduced_speed:
        for tr in st: # need to remove travel time from traces by altering starttime
            time_change = tr.stats.distance/reduced_speed
            tr.stats.starttime = tr.stats.starttime - time_change
            
        
    if endtime>starttime:
        print('Trimming')
        st.trim(starttime=starttime, endtime=endtime) 
        
    if normalize:
        st.normalize(global_max=False)
    
    rmax = np.max(r)
    l = len(st)
    max_trace_height = scale_factor/l *(km_max-km_min)
    #st.normalize(global_max=True) 
    km_per_Pa = 1.0
    last_r = -9999999
    last_asecs = -99999999
    if min_spacing_km == 0:
        min_spacing_km = 111.0  * km_max/20037.5
    max_pascals = 1500
    
    
    fh = plt.figure(figsize=figsize)
    ah = fh.add_subplot(1,1,1)
    for i, tr in enumerate(st):

        tr.detrend('linear')
        if np.max(np.abs(tr.data))>max_pascals:
            continue
        
        # distance to station and offset in km
        this_r = tr.stats.distance
        offset_km = this_r/1000
 
        # plot arrival time
        if 'arrival' in tr.stats:
            asecs = tr.stats.arrival.arrivaltime - reftime
            #if asecs < last_asecs:
            #    continue
            if plot_arrivals:
                ahours = asecs/3600
                ah.plot(ahours, offset_km, 'b.', zorder=len(st)*2+11-i)        
            last_asecs = asecs
        
        # time
        t0 = tr.stats.starttime - reftime
        t = t0 + np.arange(0, tr.stats.delta*len(tr.data), tr.stats.delta)     
        h = t/3600.0
        max_h = np.max(h) 
        
        # color
        diff_r = this_r - last_r
        brightness = 0
        if diff_r < min_spacing_km*1000:
            if diff_r < min_spacing_km/2 * 1000:
                continue
            brightness = 1.0 - (diff_r/1000)/min_spacing_km
            last_r = this_r
        else:
            last_r = this_r
        col = [brightness, brightness, brightness]
            
        # plot
        ah.plot(h, offset_km+(tr.data * km_per_Pa * scale_factor),                 color=col, linewidth=0.5, zorder=(len(st)*2+10-i)*(1.0-brightness) )

        
    if len(slopes)>0:
        print(slopes)
        for slope in slopes:
            print(slope)
            slope_speed = slope[0]
            slope_origintime = slope[1]
            if len(slope)==2:
                slope.append('r-')
            t_prime = np.array([km_min*1000/slope_speed, km_max*1000/slope_speed]) 
            d_prime = np.array([km_min, km_max])
            t_prime = d_prime*1000/slope_speed + (slope_origintime - st[0].stats.starttime)
            h_prime = t_prime/3600
            ah.plot(h_prime, d_prime, slope[2], zorder=900 )

        
    ah.set_yticks(np.arange(km_min, km_max, 1000.0))
    numhours = h[-1]-h[0]
    if numhours>20:
        stephours = np.round(numhours/20,0)
    else:
        stephours = np.round(numhours/20,1)
    #ah.set_xticks(np.arange(t[0], t[-1], 1.0))
    ah.set_xticks(np.arange(h[0], h[-1], stephours))
    ah.set_xlim(h[0],h[-1])
    ah.set_ylabel('Distance (km)')
    ah.set_xlabel('Time (hours) after %s' % reftime.strftime('%Y-%m-%d %H:%M:%S'))
    #ah.xaxis.set_major_formatter(plt.FuncFormatter('{:.0f}'.format))
    ah.yaxis.set_major_formatter(plt.FuncFormatter('{:.0f}'.format))
    #ah.autoscale(enable=True, tight=True)
    #ah.set_ylim(0, circum_earth_km/2)
    ah.set_ylim(km_min, km_max)
    ah.grid(axis='x', color = 'green', linestyle = ':', linewidth = 0.5)
    
    two_axes = True
    if two_axes:
        ah2 = ah.twinx()
        #y2 = [km_min/m_per_deg, rmax/m_per_deg]
        y2 = [1000*km_min/m_per_deg, 1000*km_max/m_per_deg]
        ah2.plot([0, 0],y2,'k-')
        ylim_km = ah.get_ylim()
        #ah2.set_ylim(ylim_km[0]*1000/m_per_deg, ylim_km[1]*1000/m_per_deg)
        ah2.set_ylabel('Distance (degrees)')
        #ah2.set_ylim(0, 180)
        ah2.set_ylim(1000*km_min/m_per_deg, 1000*km_max/m_per_deg)

    if outfile:
        fh.savefig(outfile)
    return  



def plot_combined_versus_distance(st_list, inv_list, units, corrected=False):
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(1,1,1)   
    cols = ['r','g','b']
    labels = ['IRIS LDF', 'IRIS LDO', 'RBOOM']
    for i, st_original in enumerate(st_list):
        st = st_original.copy()
        inv = inv_list[i]
        attach_station_coordinates_from_inventory(inv, st)
        attach_distance_to_stream(st, olat, olon)
        st.trim(starttime=startt+1*60, endtime=startt+18*3600) # just examine first pulse
        r = []
        a = []
        for tr in st:
            r.append(tr.stats.distance_in_km)
            #a.append(tr.stats.maxamp_corrected)
            if corrected:
                a.append(np.max(np.abs(tr.data))* np.sqrt(tr.stats.distance_in_km))
            else:
                a.append(np.max(np.abs(tr.data)))
        remove_outliers(a, f=100)     
        ax.scatter(r, a, c=cols[i], label=labels[i])
    if corrected:
        ax.set_ylabel('Reduced Pressure (%s)' % units)
    else:
        ax.set_ylabel('Maximum amplitude (%s)' % units)
    ax.set_xlabel('Distance (km)')
    ax.legend()
    plt.show()   
    return    
    
    

def remove_duplicates(st):
    st2 = obspy.core.Stream()
    for tr_original in st:
        tr = tr_original.copy()
        st2.append(tr)
        try:
            st2.merge(method=1, fill_value=0)
        except:
            st2.remove(tr)
    return st2



def machnumber(Pover, Pambient=1e5, vambient=314):
    g = 1.4 # specific heat ratio for air
    Pratio = Pover / Pambient
    f = (1 + g) * Pratio / (2 * g)
    M = np.sqrt(1 + f)
    v = vambient * M
    return (M, v)



def pressure_vs_distance(r_station, P_station, crossover_km=100, vambient=314):
    R = circum_earth_m/math.pi
    r_sphere = R * np.sin(r_station / R)
    if crossover_km == 0:
    	Preduced = P_station / (1/np.sqrt(r_sphere))
    else: 
        Preduced = P_station / (1/r_sphere + 1/np.sqrt(crossover_km * r_sphere)) 
   
    
    r = 1 # 1 km for reduced pressure
    r_list = [r]
    P_list = [Preduced]
    Pambient = 1e5
    [M,v] = machnumber(Preduced, Pambient=Pambient, vambient=vambient)
    t_list = [0]
    v_list = [v]
    M_list = [M]
    multiplier = 1.0001
    
   
    last_r = r
    while r < circum_earth_km/2: 
        r = r * multiplier
        if r > last_r + 1:
            r = last_r + 1
        r_sphere = R * np.sin(r / R)
        
        #P = Preduced/r + Preduced/np.sqrt(crossover_km * r)
        #P = Preduced/r_sphere + Preduced/np.sqrt(crossover_km * r_sphere)
        if crossover_km == 0:
    	    P = Preduced * (1/np.sqrt(r_sphere))
        else: 
            P = Preduced * (1/r_sphere + 1/np.sqrt(crossover_km * r_sphere)) 
        [M,v] = machnumber(P, Pambient=Pambient)
        r_list.append(r)
        P_list.append(P)
        v_list.append(v)
        mean_v = np.mean(v_list[-2:])
        delta_r = (r - last_r)*1000
        delta_t = delta_r/mean_v
        t_list.append(t_list[-1]+delta_t)
        M_list.append(M)
        last_r = r
    return (r_list[1:], P_list[1:], v_list[1:], t_list[1:], M_list[1:])


DATA_DIR = os.path.join(os.getenv('HOME'), 'DATA', 'HungaTonga')


####################### FUNCTIONS FOR LIGHTNING #####################

def read_wwlln():
    # Load WWLLN data
    # I got these data from Steve McNutt (email on March 1st) who got them from Bob Holzworth at UW (bobholz@uw.edu)
    # Cite: Dowden, Brundell, Rodger, Source [2002], J. Atmos. 
    # Parameters Steve asked for were -24 to -18 lat, -179 to -172 lon, 0400 to 1600 UTC
    # residual is fit residual in microseconds
    # I should add a distanceM column from -20.536, -175.382 as for the GLD360 data
    wwllnpklfile = os.path.join(DATA_DIR, 'LIGHTNING', 'WWLLN.pkl')
    if os.path.exists(wwllnpklfile):
        wwllndf = pd.read_pickle(wwllnpklfile)
    else:    
        wwllndatafile = os.path.join(DATA_DIR, 'LIGHTNING', 'A20220115.HungaTonga.loc')
        wwllndf = pd.read_csv(wwllndatafile, names=['date', 'time', 'latitude', 'longitude', 'residual', 'num_stations'])

        # Create a datetime column
        timestamp = []
        for i,row in wwllndf.iterrows():
            this_timestr = row['date'].replace('/','-')+'T'
            this_timestamp = UTCDateTime(this_timestr + row['time']).datetime
            timestamp.append(this_timestamp)
        wwllndf['datetime']=timestamp    
        wwllndf.sort_values(by='datetime',inplace=True)

        # Add other columns for plotting purposes
        wwllndf['ones']=np.ones(len(wwllndf['datetime']))
        wwllndf['cum_events']=np.cumsum(wwllndf['ones'])
        wwllndf['cum_detections']=np.cumsum(wwllndf['num_stations'])

        add_distanceM_column(wwllndf)

        # add distance weightedDetections
        weightedDetections = []
        for i, row in wwllndf.iterrows():
            try:
                weightedDetections.append(row['num_stations']/row['distanceM'])
            except:
                pass
        wwllndf['weightedDetections'] = weightedDetections       
        
        # Save to pickle files
        wwllndf.to_pickle(wwllnpklfile)        
        
    return wwllndf



def read_gld360():
    # Load Vaisala GLD360 data
    # contact is chris.vagasky@vaisala.com, see email from March 3, 2022
    # distanceM is from -20.536, -175.382. 
    gld360pklfile = os.path.join(DATA_DIR, 'LIGHTNING', 'GLD360.pkl')
    if os.path.exists(gld360pklfile):
        gld360df = pd.read_pickle(gld360pklfile)
    else:
        gld360datafile = os.path.join(DATA_DIR, 'LIGHTNING', 'Hunga_Tonga_Lightning_within_300km_December-January.csv')
        gld360df = pd.read_csv(gld360datafile)

        # Create a datetime column
        gld360df['datetime'] = [UTCDateTime(dt).datetime for dt in gld360df['time']]
        gld360df.sort_values(by='datetime',inplace=True)

        # Add other columns for plotting purposes - do not add cumulative columns unless i subset later
        gld360df['ones']=np.ones(len(gld360df['datetime']))
        gld360df['cum_events']=np.cumsum(gld360df['ones'])
        add_distanceM_column(gld360df) 
        
        # add signalPower and distance weightedPower
        signalPower = []
        weightedPower = []
        resistance = 10000 # ohms
        for i, row in gld360df.iterrows():
            try:
                thisSignalPower =  (row['signalStrengthKA']*1000)**2 * resistance
                signalPower.append(thisSignalPower)
                weightedPower.append(thisSignalPower/row['distanceM'])
            except:
                pass
        gld360df['signalPower'] = signalPower
        gld360df['weightedPower'] = weightedPower      

        # Save to pickle files
        gld360df.to_pickle(gld360pklfile)
    return gld360df



def add_distanceM_column(df):
    if not 'distanceM' in df.columns:
        df['distanceM'] = [degrees2kilometers(locations2degrees(lat, lon, sourcelat, sourcelon))*1000 for lat,lon in zip(df['latitude'],df['longitude'])]
    


# plot WWLLN data
def plot_wwlln(df, starttime='1900-01-01', endtime='2100-01-01', radiusM=12500000,                y_columns = ['cum_events', 'cum_detections']):    
    for ycol in y_columns:
        if not ycol in df.columns:
            y_columns.remove(ycol)
    df2 = df.loc[(df['datetime'] >= starttime)
            & (df['datetime'] < endtime) & df['distanceM'] < radiusM]       
    #fig, ax = plt.subplots(len(y_columns),1)
    linestyle='k.'
    
    fig = plt.figure(figsize=(12, 4))
    ax=[]
    #fig, ax = plt.subplots(len(y_columns),1)
    for i in range(len(y_columns)):
        ax.append(fig.add_subplot(len(y_columns), 1, i+1))
        #ax[i] = df2.plot.scatter(x='datetime', y=y_columns[i], c='DarkBlue', marker='.')
        ax[i].plot_date(df2['datetime'], df2[y_columns[i]], linestyle)
        ax[i].set_ylabel(y_columns[i])
        xtl = ax[i].get_xticklabels()
        plt.setp(xtl, rotation=45)
        if starttime:
            ax[i].set_xlim(starttime, endtime)
    return


        
#import matplotlib.dates as mdates 
def plot_gld360(df, starttime='1900-01-01', endtime='2100-01-01', lightning_type=None, radiusM=12500000,
                 y_columns = ['cum_events', 'signalStrengthKA', 'weightedCurrent']):
    for ycol in y_columns:
        if not ycol in df.columns:
            y_columns.remove(ycol)
    df2 = df.loc[(df['datetime'] >= starttime)
            & (df['datetime'] < endtime) & (df['distanceM'] < radiusM)] 

    linestyle = 'b.'
    if lightning_type:
        if lightning_type == 'cloud2cloud':
            df2 = df2.loc[df2['cloud']==True]
            linestyle = 'g.'
        elif lightning_type == 'cloud2ground':
            df2 = df2.loc[df2['cloud']==False]
            linestyle = 'r.'
            
    df2['cum_events'] = np.cumsum(df2['ones'])
    
    st = obspy.Stream()
    fig = plt.figure(figsize=(12, 4))
    ax=[]
    #fig, ax = plt.subplots(len(y_columns),1)
    for i in range(len(y_columns)):
        ax.append(fig.add_subplot(len(y_columns), 1, i+1))
        #ax[i] = df2.plot.scatter(x='datetime', y=y_columns[i], c='DarkBlue', marker='.')
        ax[i].plot_date(df2['datetime'], df2[y_columns[i]], linestyle)
        ax[i].set_ylabel(y_columns[i])
        xtl = ax[i].get_xticklabels()
        plt.setp(xtl, rotation=45)
        #ax[i].xaxis.set_major_locator(mdates.DayLocator(interval=4))
        #ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%d'))
        if starttime:
            ax[i].set_xlim(starttime, endtime)  
        #tr = obspy.Trace() 
        #tr.data = 
    return



def filter_df(df, starttime, endtime, radiusM=12500000 ): # filter a lightning data frame by time or distance
    df2 = df.copy()
    
    # time filter
    df3 = df2.loc[(df2['datetime'] >= starttime) & (df2['datetime'] < endtime)]
    
    # spatial filter
    df4 = df3.loc[ (df3['distanceM'] < radiusM)] 
    
    return df4



# Next we want to bin the data in 0.01 s intervals.
def bin_by_seconds(dfin, dt=60):    
    df2 = dfin.set_index('datetime', inplace=False)
    counts = df2.resample('%ds' % dt).sum()
    #print(counts)
    return counts



def column2Trace(df, colname, seconds, net='', sta='', chan=''):
    tr = obspy.Trace()
    tr.id = '%s.%s.%d.%s' % (net, sta, seconds, chan)
    print(df)
    tr.stats.starttime = obspy.UTCDateTime(df.index[0])
    tr.stats.delta = seconds
    tr.data = np.array(df[colname])
    return tr
