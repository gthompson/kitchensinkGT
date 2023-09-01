import os
import numpy as np
from obspy import Stream
import obspyGT.SDS
from obspyGT.InventoryTools import attach_station_coordinates_from_inventory, attach_distance_to_stream


class RSAMobj:
    
    def __init__(self, st=None, inv=None, sampling_interval=60.0, verbose=False, metric='mean', freqmin=0.1, zerophase=False, corners=2 ):
        self.stream = Stream()
        self.metric = metric
        if isinstance(st, Stream):
            st.detrend(type='linear') # remove a linear trend
            # filter >0.5 Hz by default for RSAM
            if inv:
                st.filter('highpass', freq=freqmin, zerophase=False, corners=corners) 
                st.remove_response(output='VEL', inventory=inv)    
            else:
                print('No inventory for RSAM calculation. Only detrended.')  
            for tr in st:
                this_tr = tr.copy()
                x = tr.data # get the data
                # now we want to reshape the data vector into an array, so we can take advantage of np.mean()
                s = np.size(x) # find the size of the data vector
                nc = int(tr.stats.sampling_rate * sampling_interval) # number of columns
                nr = int(s / nc) # number of rows
                x = x[0:nr*nc] # cut off any trailing samples
                y = x.reshape((nr, nc))
                if verbose:
                    print('%s: size %d' % (tr.id, s))
                    print('%s: reshaped to %d x %d array (%d samples)' % (tr.id, nr, nc, nr * nc))
                if metric=='mean':
                    this_tr.data = np.nanmean(abs(y),axis=1) # compute mean for each row (this is vectorised; faster than a for loop)
                if metric=='median':
                    this_tr.data = np.nanmedian(abs(y),axis=1) # compute mean for each row (this is vectorised; faster than a for loop)
                if metric=='max':
                    this_tr.data = np.nanmax(abs(y),axis=1) # compute mean for each row (this is vectorised; faster than a for loop)
                if metric=='rms':
                    this_tr.data = np.rms(abs(y),axis=1) # compute mean for each row (this is vectorised; faster than a for loop)     
                this_tr.stats.sampling_rate = 1.0 / sampling_interval # update the sampling rate
                self.stream.append(this_tr)
        self.sampling_interval = sampling_interval

    def plot(self, equal_scale=False, type='normal', percentile=None):
        if type=='normal':
            normalplot(self.stream, equal_scale=equal_scale, percentile=percentile)


    def _fix_traceid(self):
        # See: https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/
        for tr in self.stream:
            if tr.stats.delta == 1:
                tr.stats.channel = 'L' + tr.stats.channel[1:]
            if tr.stats.delta == 60:
                tr.stats.channel = 'V' + tr.stats.channel[1:]
            if tr.stats.delta == 600:
                tr.stats.channel = 'U' + tr.stats.channel[1:]          


    def write(self, SDS_TOP):
        RSAMSDS_TOP = os.path.join(SDS_TOP,'RSAM',self.metric)
        #self._fix_traceid()
        RSAMSDSobj = obspyGT.SDS.SDSobj(RSAMSDS_TOP, streamobj=self.stream)
        RSAMSDSobj.write()           
        
    
    def read(self, startt, endt, SDS_TOP, metric='mean', speed=2):
        RSAMSDS_TOP = os.path.join(SDS_TOP,'RSAM',self.metric)
        thisSDSobj = obspyGT.SDS.SDSobj(RSAMSDS_TOP)
        thisSDSobj.read(startt, endt, speed=speed)
        print(thisSDSobj.stream)
        self.stream = thisSDSobj.stream
        self.metric=metric
        self.sampling_interval=self.stream[0].stats.delta



class ReducedDisplacementObj():
    # Difference from RSAM is 
    # (1) Correct to displacement first 
    # (2) Filter, 
    # (3) Apply geometric spreading 
    # (4) Correct for Q
    def __init__(self, st=None, inv=None, sampling_interval=60.0, verbose=False, metric='median', \
            freqmin=0.1, freqmax=15.0, zerophase=False, corners=2, peakf=2.0, wavespeed=2000, centerlat=None, centerlon=None ):
        self.stream = Stream()
        self.metric = metric
        if isinstance(st, Stream):
            st.detrend(type='linear') # remove a linear trend
            # filter 0.5-15.0 by default for Drs
            if inv:
                st.filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=False, corners=corners) 
                st.remove_response(output='DISP', inventory=inv) 
            else:
                print('No inventory for DRS calculation. Cannot compute.') 
                return 
            attach_station_coordinates_from_inventory(inv, st)
            attach_distance_to_stream(st, centerlat, centerlon)
            for tr in st:
                this_tr = tr.copy()
                x = (tr.data*100) * np.sqrt(tr.stats.distance*100 * peakf * wavespeed * 100) # everything in cm
                # now we want to reshape the data vector into an array, so we can take advantage of np.mean()
                s = np.size(x) # find the size of the data vector
                nc = int(tr.stats.sampling_rate * sampling_interval) # number of columns
                nr = int(s / nc) # number of rows
                x = x[0:nr*nc] # cut off any trailing samples
                y = x.reshape((nr, nc))
                if verbose:
                    print('%s: size %d' % (tr.id, s))
                    print('%s: reshaped to %d x %d array (%d samples)' % (tr.id, nr, nc, nr * nc))
                if metric=='mean':
                    this_tr.data = np.nanmean(abs(y),axis=1) # compute mean for each row (this is vectorised; faster than a for loop)
                if metric=='median':
                    this_tr.data = np.nanmedian(abs(y),axis=1) # compute mean for each row (this is vectorised; faster than a for loop)
                if metric=='max':
                    this_tr.data = np.nanmax(abs(y),axis=1) # compute mean for each row (this is vectorised; faster than a for loop)
                if metric=='rms':
                    this_tr.data = np.rms(abs(y),axis=1) # compute mean for each row (this is vectorised; faster than a for loop)     
                this_tr.stats.sampling_rate = 1.0 / sampling_interval # update the sampling rate
                self.stream.append(this_tr)
        self.sampling_interval = sampling_interval

    def plot(self, equal_scale=False, type='normal', percentile=None):
        if type=='normal':
            normalplot(self.stream, equal_scale=equal_scale, percentile=percentile) 

        elif type=='iceweb':
            import matplotlib.pyplot as plt
            import matplotlib.dates as mdates
            from math import log2
            plt.rcParams["figure.figsize"] = (10,6)
            fig, ax = plt.subplots()
            for tr in self.stream:
                x = tr.data
                 # now we want to reshape the data vector into an array, so we can take advantage of np.mean()
                s = np.size(x) # find the size of the data vector
                nc = np.max((1, int(log2(s/260))))  # number of columns
                nc = nc * 4
                nr = int(s / nc) # number of rows
                x = x[0:nr*nc] # cut off any trailing samples
                y = x.reshape((nr, nc))
                y2 = np.nanmax(y,axis=1)    
                t = tr.times("utcdatetime")[::nc]   
                t = [this_t.datetime for this_t in t] 
                #print(t) 
                t = t[:len(y2)]      
                ax.semilogy(t, y2,'.', label='%s' % tr.id) #, alpha=0.03)
            ax.format_xdata = mdates.DateFormatter('%H')
            ax.legend()
            plt.xticks(rotation=45)
            plt.ylim((0.2, 100)) # IceWeb plots went from 0.05-30
            plt.yticks([0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0], \
                ['0.2', '0.5', '1', '2', '5', '10', '20', '50', '100'])
            plt.ylabel(r'$D_{RS}$ ($cm^{2}$)')
            plt.xlabel(r'UTC / each point is max $D_{RS}$ in %d minute window' % (tr.stats.delta * nc /60))
            plt.title('Reduced Displacement (%s)\n%s to %s' % (r'$D_{RS}$', t[0].strftime('%d-%b-%Y %H:%M:%S UTC'), t[-1].strftime('%d-%b-%Y %H:%M:%S UTC')))
        return 0


    def _fix_traceid(self):
        # See: https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/
        for tr in self.stream:
            if tr.stats.delta == 1:
                tr.stats.channel = 'L' + tr.stats.channel[1:]
            if tr.stats.delta == 60:
                tr.stats.channel = 'V' + tr.stats.channel[1:]
            if tr.stats.delta == 600:
                tr.stats.channel = 'U' + tr.stats.channel[1:]          


    def write(self, SDS_TOP):
        DrsSDS_TOP = os.path.join(SDS_TOP,'DRS',self.metric)
        #self._fix_traceid()
        DrsSDSobj = obspyGT.SDS.SDSobj(DrsSDS_TOP, streamobj=self.stream)
        DrsSDSobj.write()                
    
    def read(self, startt, endt, SDS_TOP, metric='mean', speed=2):
        DrsSDS_TOP = os.path.join(SDS_TOP,'DRS',self.metric)
        thisSDSobj = obspyGT.SDS.SDSobj(DrsSDS_TOP)
        thisSDSobj.read(startt, endt, speed=speed)
        print(thisSDSobj.stream)
        self.stream = thisSDSobj.stream
        self.metric=metric
        self.sampling_interval=self.stream[0].stats.delta  

def daily_wrapper(startt, endt, SDS_TOP, centerlat, centerlon, searchRadiusDeg, \
        fdsnURL="http://service.iris.edu", channel='BHZ', freqmin=0.5, freqmax=15.0, \
        zerophase=False, corners=2, sampling_interval=60.0, writeSDS=True, writeDRS=True):
    import obspyGT.SDS
    import obspyGT.FDSNtools
    import obspyGT.InventoryTools
    secsPerDay = 86400  
    while startt<endt:
        print(startt)
        eod = startt+secsPerDay-1/10000
        # read from SDS - if no data download from FDSN

        thisSDSobj = obspyGT.SDS.SDSobj(SDS_TOP) 
        inv = obspyGT.FDSNtools.get_inventory(fdsnURL, startt, eod, centerlat, centerlon, \
            searchRadiusDeg, channel=channel, overwrite=False )
        if thisSDSobj.read(startt, eod, speed=2): # non-zero return value means no data in SDS so we will use FDSN
            # read from FDSN
            trace_ids = obspyGT.InventoryTools.inventory2traceid(inv)
            st = obspyGT.FDSNtools.get_stream(fdsnURL, trace_ids, startt, eod, overwrite=True)
            thisSDSobj.stream = st
            if writeSDS:
                thisSDSobj.write() # save raw data to SDS

        # compute instrument-corrected RSAM
        thisRSAMobj = RSAMobj(st=thisSDSobj.stream.copy(), inv=inv, sampling_interval=sampling_interval, \
                              freqmin=freqmin, zerophase=zerophase, corners=corners)
        thisRSAMobj.write(SDS_TOP) # write RSAM to an SDS-like structure
        
        # compute/write reduced displacement
        if writeDRS:
            thisDRSobj = ReducedDisplacementObj(st=thisSDSobj.stream.copy(), inv=inv, sampling_interval=sampling_interval, \
                                freqmin=freqmin, freqmax=freqmax, zerophase=zerophase, corners=corners, \
                                     centerlat=centerlat, centerlon=centerlon )
            thisDRSobj.write(SDS_TOP) # write Drs to an SDS-like structure
    

        startt+=secsPerDay # add 1 day 

def normalplot(st, equal_scale=False, percentile=None):
    hf = st.plot(handle=True, equal_scale=equal_scale); # standard ObsPy plot
    # change the y-axis so it starts at 0
    allAxes = hf.get_axes()
    ylimupper = [ax.get_ylim()[1] for ax in allAxes]
    print(ylimupper)
    if percentile:
        ylimupper = np.array([np.percentile(tr.data, percentile) for tr in st])*1.1
    # if equal_scale True, we set to maximum scale
    print(ylimupper)
    ymax=max(ylimupper)
    for i, ax in enumerate(allAxes):
        if equal_scale==True:
            ax.set_ylim([0, ymax])
        else:
            ax.set_ylim([0, ylimupper[i]])  

if __name__ == "__main__":
    print("This is the RSAM module. It also handles Reduced Displacement")
    print('Example: RSAM_wrapper_Shishaldin.ipynb on 2023/08/31')