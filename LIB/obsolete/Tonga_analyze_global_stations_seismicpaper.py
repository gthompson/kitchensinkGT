#!/usr/bin/env python
# coding: utf-8
# Functions from Tonga_analyze_global_stations_seismicpaper.ipynb
#import os, sys, obspy
#import numpy as np
#get_ipython().run_line_magic('matplotlib', 'inline')
#import matplotlib.pyplot as plt
#import pickle
LIBpath = os.path.join( os.getenv('HOME'),'src','kitchensinkGT', 'LIB')
sys.path.append(LIBpath)
from libseisGT import inventory2traceid, get_FDSN_inventory, clean_trace, \
    attach_station_coordinates_from_inventory
from obspy.clients.fdsn import Client
from obspy.geodetics import locations2degrees, degrees2kilometers, kilometers2degrees
#from scipy.stats import linregress

def get_fdsn_identifier(fdsnURL):
    prepend=''
    if 'iris' in fdsnURL:
        prepend='iris_'
    if 'shake' in fdsnURL:
        prepend='rboom_'
    if 'geonet' in fdsnURL:
        prepend='nz_'  
    return prepend

def get_inventory_and_waveforms(fdsnClient, searchRadiusDeg, olat, olon, startt, endt, chanstring, data_root, overwrite=False, network='*', station='*', location='*'):
    fdsnClient = Client(fdsnURL)  
    if not os.path.isdir(data_root):
        print('No such directory %s' % data_root)
        return None, None
    
    # STEP 1 - FIND STATIONS BY REQUESTING AN INVENTORY
    prepend = get_fdsn_identifier(fdsnURL)
    stationXmlFile = os.path.join(data_root, '%s%s_within_%.0f_degrees.xml' % (prepend, chanstring, searchRadiusDeg))
    mseedfile = stationXmlFile.replace('.xml','.mseed')
    print(stationXmlFile, mseedfile)
    
    if os.path.exists(stationXmlFile) and overwrite==False:
        inv = obspy.read_inventory(stationXmlFile)
    else:
        inv = get_FDSN_inventory(fdsnClient, startt, stationXmlFile, network, olat, olon, \
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
            st.write(stationXmlFile.replace('.xml','.mseed'))
        
    return st, inv

def reconstitute_stream(st, inv, fmin=0.0001):
    T=1.0/fmin
    
    # deal with multiple traces with same trace ID, by appending and merging one at a time,
    # and removing any time there is a fail.
    st2 = obspy.core.Stream()
    for tr_original in st:
        tr = tr_original.copy()
        st2 = st2.append(tr)
        try:
            st2.merge(method=1, fill_value=0)
        except:
            st2.remove(tr)
            pass
        
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
                tr.remove_response(inventory=inv, pre_filt=pre_filt, output='DISP')
            else:
                tr.remove_response(inventory=inv, pre_filt=pre_filt )
            reconstituted.append(tr)
            #print(tr.id, ' reconstituted')
        except:
            #print(tr.id, ' NOT reconstituted')
            failed_ids.append(tr.id)
            pass
        tr.detrend('linear')
        tr.stats['maxamp_corrected'] = np.max(np.abs(tr.data))
    print('Failed to reconstitute: ', failed_ids)

    #reconstituted.plot(equal_scale=True);
    
    return reconstituted


def attach_station_coordinates_from_inventory(inventory, st):
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
    
def plot_amplitude_versus_distance(st, units, reftime=None, cmin=300, cmax=420, duration=7200):
    r = []
    a = []
    t = []
    for tr in st:
        
        if reftime:
            this_r = tr.stats.distance
            Pmax = []
            tmax = []
            stime = reftime + this_r/cmax
            etime = reftime + this_r/cmin + duration
            if etime < tr.stats.endtime:
                tr0=tr.copy().trim(starttime=stime, endtime=endt)
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
                tr0=tr.copy().trim(starttime=stime, endtime=endt)
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
    radius_earth = circum_earth_m/(2*3.1415)
    for this_r in r2: 
        real_r = radius_earth * np.arcsin(this_r/radius_earth)
        print(this_r, real_r)
        r3.append(real_r)
        y2.append(PR / np.power(real_r, 0.5))
    y2 = np.array(y2)
    r2 = np.array(r2)
    r3 = np.array(r3)
    #print(y2)
    
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(1,1,1)
    ax.loglog(r, a, 'o')
    ax.loglog(r2, y, 'r:')
    ax.loglog(r2, y2, 'g:')
    ax.set_ylim(10,100000)
    ax.set_xlim(1,42000)
    ax.set_ylabel('Maximum amplitude (%s)' % units)
    ax.set_xlabel('Distance (km)')
    plt.grid(which='both')
    #plt.savefig('reduced_pressure.eps')
    
    fig = plt.figure()
    plt.plot(r2, r3/r2, 'o')
    
    if False:
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(1,1,1)
        ax.scatter(t, r)
        ax.set_xlabel('Maximum time (s)')
        ax.set_ylabel('Distance (km)')
        plt.show()    
    
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
            
            # only recommend to keep a trace if amplitude less than twice previous
            this_max_amp = np.max(np.abs(tr.data))
            print('max amp = ', this_max_amp)
            this_distance = tr.stats.distance
            if this_distance > rmax: # skip if greater than rmax
                continue
            print('distance = ', this_distance)
            recommend = 'n'
            if this_max_amp < 2 * last_max_amp:
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

def analyze_clientchan(fdsnClient, chanstring, fmin=0.001, network='*', station='*', location='*'):
    T = np.round(1/fmin,0)
    
    this_stream, this_inventory = get_inventory_and_waveforms(fdsnClient, searchRadiusDeg, olat, olon, startt, endt, chanstring, DATA_ROOT, network=network, station=station, location=location, overwrite=False)
    attach_station_coordinates_from_inventory(this_inventory, this_stream)
    attach_distance_to_stream(this_stream, olat, olon)
    
    prepend = get_fdsn_identifier(fdsnClient)
    this_reconstituted = reconstitute_stream(this_stream, this_inventory, fmin=fmin)
    reconfile = os.path.join(DATA_ROOT, '%s%s_reconstituted_%.0f.mseed' % (prepend, chanstring, T))    
    #if chanstring[0]!='L':
    #    median_downsample_to_1Hz(this_reconstituted)
    #    reconfile = os.path.join(DATA_ROOT, '%s%s_reconstituted_%.0f_downsampled.mseed' % (prepend, chanstring, T))
    if not os.path.exists(reconfile):
        this_reconstituted.write(reconfile)
    
    good_ids_pickle = os.path.join(DATA_ROOT, '%sgood_%s_ids_within_%.0f_degrees.pkl' % (prepend, chanstring, searchRadiusDeg))
    if os.path.exists(good_ids_pickle):
        with open(good_ids_pickle, 'rb') as f:
            good_ids = pickle.load(f)        
        this_good = select_previously_manually_selected_traces(this_reconstituted, good_ids)
    else:
        this_good = manually_select_good_traces(this_reconstituted)
        good_ids = []
        for tr in this_good:
            good_ids.append(tr.id)
        print(good_ids)
        print(len(good_ids))
        with open(good_ids_pickle, 'wb') as f:
            pickle.dump(good_ids, f)    
        this_good.write(os.path.join(DATA_ROOT, '%s%s_good.mseed' % (prepend, chanstring)))   
    
    if chanstring[1]=='D':
        units = 'Pa'
        plot_amplitude_versus_distance(this_good, units)
        plot_reduced_pressure(this_good, units)
    #elif chanstring[1]=='H':
    #    units = 'm/s'
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

def get_distance_vector(st):
    r = []
    for tr in st:
        r.append(tr.stats.distance)
    return r
        
    
def order_traces_by_distance(st, r): 
    st2 = obspy.core.Stream()
    indices = np.argsort(r)
    #print(indices)
    for i in indices:
        tr = st[i].copy()
        st2.append(tr)
        #print(tr)
    return st2

def pick_pressure_onset(st, otime, cmin, cmax):
    print('Pick onset times with left mouse button. Middle button to skip trace.')
    for tr in st:
        mintime = otime + tr.stats.distance/cmax
        maxtime = otime + tr.stats.distance/cmin
        winsecs = maxtime - mintime
        
        good_trace = False
        
        while winsecs > 180:
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
                    print(atime)
                    winsecs = winsecs / 5
                    mintime = atime - winsecs/2
                    maxtime = atime + winsecs/2 
                else:
                    break
            else:
                break
        
        
        if good_trace:
            traveltime = atime - otime
            tr.stats['arrival'] = obspy.core.util.AttribDict({
                                'arrivaltime':atime,
                                'c':tr.stats.distance/traveltime,
                                'traveltime':traveltime})
            print(tr.stats.arrival)
        

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

def plot_record_section(st2, reftime, outfile=None, starttime=0, endtime=0, \
                        km_min=0, km_max=40075.0/2, slopes = [], \
                        scale_factor=1.0, plot_arrivals=False, \
                        figsize=(16,16), min_spacing_km=0, reduced_speed=None, normalize=False ):
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
        ah.plot(h, offset_km+(tr.data * km_per_Pa * scale_factor), \
                color=col, linewidth=0.5, zorder=(len(st)*2+10-i)*(1.0-brightness) )

        
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

def process_1(st, otime=otime, cmin=300, cmax=430, decimation_factor=10):
    st_downsample = st.copy()
    st_downsample.decimate(decimation_factor)
    r = get_distance_vector(st_downsample)
    print(r)
    st_downsample_pick = order_traces_by_distance(st_downsample, r) 
    st_downsample_pick = manually_select_good_traces(st_downsample_pick, trim_mins=0, record_section=False, auto=False, wlen=7200)
    #pick_pressure_onset(st_downsample_pick, otime, cmin, cmax)
    return st_downsample_pick
