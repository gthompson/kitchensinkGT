import numpy as np
import matplotlib.pyplot as plt
import obspy.core
from scipy import stats
from obspy.core.stream import Stream
from obspy.core import UTCDateTime
import types

def hanningcosine(data, fraction=1/3):
    """ create a cosine tapered hanning window with a flat central section """
    from numpy import round
    numsamplesdata = len(data)
    numsamplestaper = fraction * round(numsamplesdata)
    from scipy import signal
    h = signal.hanning(numsamplestaper * 2)
    h1 = h[:numsamplestaper]
    h2 = h[-numsamplestaper]
    win = [h1, [1.0]*numsamplesdata, h2]
    return win

def nanrms(x, axis=None):
    return np.sqrt(stats.nanmean(x**2, axis=axis))

def patch(self):

    def summary(st):
        print "STAT.CHA |  max      rms      std      offset "
        print "---------|------------------------------------"
        for i in range(len(st)):
             # We could detrend, but in case of spikes, subtracting the median may be better
             #st[i].detrend()
             offset = np.median(st[i].data)
             y = st[i].data - offset
             ymax = np.max(np.abs(y))
             yrms = nanrms(y)
             ystd = stats.nanstd(y)
             print "%s.%s | %8.1e %8.1e %8.1e %8.1e" % (st[i].stats.station, st[i].stats.channel, ymax, yrms, ystd, offset)
    st.summary = types.MethodType(summary, st)

    def mulplt(st, bottomlabel='', ylabels=[]):
        plt.figure()
        # start time as a Unix epoch
        startepoch = st[0].stats.starttime.timestamp
        # create empty set of subplot handles - without this any change to one affects all
        axh = []
        # loop over all stream objects
        n = len(st)
        for i in range(n):
            # add new axes handle for new subplot
            #axh.append(plt.subplot(n, 1, i+1, sharex=ax))
            axh.append(plt.subplot(n, 1, i+1))
            # time vector, t, in seconds since start of record section
            t = np.linspace(st[i].stats.starttime.timestamp - startepoch,
                st[i].stats.endtime.timestamp - startepoch,
                st[i].stats.npts)
            # We could detrend, but in case of spikes, subtracting the median may be better
            #st[i].detrend()
            offset = np.median(st[i].data)
            y = st[i].data - offset
            # PLOT THE DATA
            axh[i].plot(t, y)
            # remove yticks because we will add text showing max and offset values
            axh[i].yaxis.set_ticks([])
            # remove xticklabels for all but the bottom subplot
            if i < n-1:
                axh[i].xaxis.set_ticklabels([])
            else:
                # for the bottom subplot, also add an xlabel with start time
                if bottomlabel=='':
                    plt.xlabel("Starting at %s" % (st[0].stats.starttime) )
                else:
                    plt.xlabel(bottomlabel)
            # default ylabel is station.channel
            if ylabels==[]:
                plt.ylabel(st[i].stats.station + "." + st[i].stats.channel, rotation=0)
            else:
                plt.ylabel(ylabels[i])
            # explicitly give the maximum amplitude and offset(median)
            plt.text(0, 1, "max=%.1e offset=%.1e" % (np.max(np.abs(y)), offset),
                horizontalalignment='left',
                verticalalignment='top',transform=axh[i].transAxes)
        # change all font sizes
        plt.rcParams.update({'font.size': 8})
        # show the figure
        plt.show()
    st.mulplt = types.MethodType(mulplt,st)

    def to_sam(st, metric='mean'):
        """
        take waveform data and compute 1 minute RSAM (or a similar metric)
        Should I just be decimating or resampling the data?
        metric can be 'mean', 'median', 'std', 'rms', 'max', 'energy'
        """
        from sam import sam
        samobject = list()
        for tr in st:
            data = list()
            wdataset = eend = None
            # create the data variable 1 minute at a time
            # SCAFFOLD: I should find minute marks using floor(starttime) and ceil(endtime)
            numminutes = int((len(tr.data)/tr.stats.sampling_rate)/60)
            samples_per_minute = int(tr.stats.sampling_rate * 60)
            for i in range(numminutes):
                pstart = i * samples_per_minute
                pend = (i+1) * samples_per_minute - 1
                wdatasubset = tr.data[pstart:pend]
                cmd = 'data.append(stats.nan' + metric + '(wdatasubset))'
                eval(cmd)
            #cmd2 = 'stats.nan' + metric + '(data)'
            #print eval(cmd2)
            # create the samobject for this trace
            samobject.append(sam(metric, tr.stats.starttime, tr.stats.endtime, 
                tr.stats.station, tr.stats.channel, data))
        return samobject
    st.to_sam = types.MethodType(to_sam, st)

    ######################################################################
    ##### REWRITE THESE IN TERMS OF OBSPY STREAM METHODS #################
    ##### ALSO THINK ABOUT GAPS, SPIKES, VLPS            #################
    ######################################################################
    def reconstitute(self, deconvolve_flag=False, integrate_flag = False):
        """
        reconstitute waveform
        first filter it robustly between 0.8 and 15 Hz
        then apply calibration correction (or deconvolve if deconvolve_flag)
        """
        self.detrend()
        self.filter(0.8, 15, 4)
        if deconvolve_flag:
                self.deconvolve()
        else:
                self.calibrate()
        if integrate_flag:
                self.integrate()
                self.seismogramtype = 'displacement'
        else:
                self.seismogramtype = 'velocity'

    def filter(self, lowcut=0.8, highcut=15.0, poles=4):
        """
        filter
        # create d = [wdata-reversed][wdata][wdata-reversed]
        # apply 50% hanning window
        # apply strong bandpass filter 0.8 - 15 Hz, 2 way
        # wdata = extract middle 3rd
        """
        wdata = self.data
        numsamples = len(wdata)
        wdata3 = [fliplr(wdata), wdata, fliplr(wdata)]
        win = hanningcosine(wdata3)
        wdata3 = wdata3 * win
        # create bandpass filter lowcut to highcut
        # apply both ways on wdata3
        # extract middle 3rd
        wdata = wdata3[numsamples + 1: numsamples * 2]
        return wdata

if __name__ == "__main__":

    from obspy.fdsn import Client
    client = Client()
    t = obspy.core.UTCDateTime("2010-02-27T06:45:00.000")
    st = client.get_waveforms("IU", "ANMO", "00", "BHZ", t, t + 60 * 5)
    print st[0]
    print st[0].stats
    patch(st)
    #st.mulplt()
    st.summary()
    st.detrend()
    st.summary()
    st.plot()
    s = st.to_sam()
    for i in s:
        print i
        i.plot()
    for i in s:
        print i
        i.plot()
