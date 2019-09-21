import numpy as np
import matplotlib.pyplot as plt
import obspy.core
from scipy import stats
from obspy.core.stream import Stream
from __future__ import print_function
from __future__ import unicode_literals
from future import standard_library  # NOQA
from future.builtins import zip
from future.builtins import range
from future.builtins import open
from future.builtins import str
from future.utils import native_str
from glob import glob, has_magic
from obspy.core.trace import Trace
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.util import NamedTemporaryFile
from obspy.core.util.decorator import map_example_filename
from obspy.core.util.base import ENTRY_POINTS, _readFromPlugin, \
    _getFunctionFromEntryPoint
from obspy.core.util.decorator import uncompressFile, raiseIfMasked
from obspy.core import compatibility
from pkg_resources import load_entry_point
import pickle
import copy
import fnmatch
import math
import numpy as np
import os
import warnings


@map_example_filename("pathname_or_url")
def read(pathname_or_url=None, format=None, headonly=False, starttime=None,
         endtime=None, nearest_sample=True, dtype=None, apply_calib=False,
         **kwargs):
    # add default parameters to kwargs so sub-modules may handle them
    kwargs['starttime'] = starttime
    kwargs['endtime'] = endtime
    kwargs['nearest_sample'] = nearest_sample
    # create stream
    st = Stream()
    if pathname_or_url is None:
        # if no pathname or URL specified, return example stream
        st = _createExampleStream(headonly=headonly)
    elif not isinstance(pathname_or_url, (str, native_str)):
        # not a string - we assume a file-like object
        pathname_or_url.seek(0)
        try:
            # first try reading directly
            stream = _read(pathname_or_url, format, headonly, **kwargs)
            st.extend(stream.traces)
        except TypeError:
            # if this fails, create a temporary file which is read directly
            # from the file system
            pathname_or_url.seek(0)
            with NamedTemporaryFile() as fh:
                fh.write(pathname_or_url.read())
                st.extend(_read(fh.name, format, headonly, **kwargs).traces)
        pathname_or_url.seek(0)
    elif "://" in pathname_or_url:
        # some URL
        # extract extension if any
        suffix = os.path.basename(pathname_or_url).partition('.')[2] or '.tmp'
        with NamedTemporaryFile(suffix=suffix) as fh:
            fh.write(compatibility.urlopen(pathname_or_url).read())
            st.extend(_read(fh.name, format, headonly, **kwargs).traces)
    else:
        # some file name
        pathname = pathname_or_url
        for file in sorted(glob(pathname)):
            st.extend(_read(file, format, headonly, **kwargs).traces)
        if len(st) == 0:
            # try to give more specific information why the stream is empty
            if has_magic(pathname) and not glob(pathname):
                raise Exception("No file matching file pattern: %s" % pathname)
            elif not has_magic(pathname) and not os.path.isfile(pathname):
                raise IOError(2, "No such file or directory", pathname)
            # Only raise error if no starttime/endtime has been set. This
            # will return an empty stream if the user chose a time window with
            # no data in it.
            # XXX: Might cause problems if the data is faulty and the user
            # set starttime/endtime. Not sure what to do in this case.
            elif not starttime and not endtime:
                raise Exception("Cannot open file/files: %s" % pathname)
    # Trim if times are given.
    if headonly and (starttime or endtime or dtype):
        msg = "Keyword headonly cannot be combined with starttime, endtime" + \
            " or dtype."
        warnings.warn(msg, UserWarning)
        return st
    if starttime:
        st._ltrim(starttime, nearest_sample=nearest_sample)
    if endtime:
        st._rtrim(endtime, nearest_sample=nearest_sample)
    # convert to dtype if given
    if dtype:
        for tr in st:
            tr.data = np.require(tr.data, dtype)
    # applies calibration factor
    if apply_calib:
        for tr in st:
            tr.data = tr.data * tr.stats.calib
    return st

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

class StreamGT(Stream):
    def __init__(self):
        Stream.__init__(self)
    
    #def read(self, ):
    #    Stream.read(self)

    def statistics(self):
        #print "STAT.CHA |  max      rms      std      offset "
        #print "---------|------------------------------------"
        for i in range(len(st)):
             # We could detrend, but in case of spikes, subtracting the median may be better
             #st[i].detrend()
             offset = np.median(st[i].data)
             y = st[i].data - offset
             ymax = np.max(np.abs(y))
             yrms = nanrms(y)
             ystd = stats.nanstd(y)
             print "%s.%s | %8.1e %8.1e %8.1e %8.1e" % (st[i].stats.station, st[i].stats.channel, ymax, yrms, ystd, offset)

    def nanrms(x, axis=None):
        return np.sqrt(stats.nanmean(x**2, axis=axis))

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

    def to_sam(st, method='rsam'):
        """
        take waveform data and compute 1 minute RSAM (or a similar metric)
        Should I just be decimating or resampling the data?
        """
        import sam
        samobject = list()
        for tr in st:
            data = list()
            wdataset = eend = None
            # create the data variable 1 minute at a time
            # SCAFFOLD: I should find minute marks using floor(starttime) and ceil(endtime)
            numminutes = (len(tr.data)/tr.samprate)/60
            samples_per_minute = tr.samprate * 60
            for i in range(numminutes):
                pstart = i * tr.samprate
                pend = (i+1) * tr.samprate - 1
                wdatasubset = tr.data[pstart:pend]
                if method == 'esam':
                    data[i] = nansum(wdatasubset)
                else:
                    data[i] = nanmedian(wdatasubset)
            # create the samobject for this trace
            eend = tr.estart + len(tr.data)/tr.samprate # I think there is a stream method
                                                            # that returns this
            samobject.append(sam.sam(method, tr.estart, eend, tr.sta, tr.chan, '', data))
        return samobject

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
    st = client.get_waveforms("IU", "ANMO", "00", "BHZ", t, t + 60 * 60)
    #st.plot() 
    st.write('iu.anmo.00.bhz', 'sac')
    
    st2 = StreamGT()
    st2.read('iu.anmo.00.bhz', 'sac')
    st2.mulplt()
