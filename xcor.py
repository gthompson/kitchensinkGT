""" xcor.py

This script contains some simple methods to cross-correlate template waveforms using the 
normalized Pearson Cross-correlation algorithm.  It is designed to work with Obspy (obspy.org)
time-series data, but is generic enough to work with any time-series data.  It works well
for event and phase identification.  While cross-correlation is not new and I do not have a
substantive paper on the matter.  I did spend a lot of time vectorizing the algorithm and
would certainly appreciate recognition for this effort.  If you find this code helpful, I
would appreciate it if you cited my paper for which this was developed.  I have omitted
the slow version of cross correlation since I can't imagine anyone would want to work 
slowly.  The vectorized version does have a drawback you are limited on the number of points
your time-series can contain before you use too much memory.  The memory usage is almost
instantaneous and gains amazing computational performance.  Generally this time is more 
than an hour of time-series data at 200 Hz, but may vary.

Dependencies:
obspy.core
obspy.signal
numpy
matplotlib

Example of usage:
# This is a simple example but can be run without any data in hand
from obspy.core import *
from obspy.signal import *
from xcor import *


# Variables to control the behavior
bandpass=[1.0,5.0]  # Define our bandpass min and max values
taper=0.1         # Percent taper to apply to the template
tmplt_dur=12.0 # Duration to use for template from origintime in seconds
xcor_cut=0.7     # Cross correlation value sufficient to identify

# I generally build my templates starting at the origintime of well identified
# events and then make sure I capture both the P and S wave
# Phase correlations can be a bit more challenging depending on the data with which you 
# are working. 
st=read()
tw=st[0]
# where tmplt.time is an origintime I use a sqlalchemy class to manage data
tw.trim(starttime=UTCDateTime(2009,8,24,0,20,6),endtime=UTCDateTime(2009,8,24,0,20,3)+tmplt_dur)
tw.filter("bandpass",freqmin=bandpass[0],freqmax=bandpass[1],zerophase=True)
tw.detrend('constant')
tap=cosTaper(tw.stats.npts,taper)
tw.data=tw.data*tap
st=read()
tr=st[0]
tr.filter("bandpass",freqmin=bandpass[0],freqmax=bandpass[1],zerophase=True)
if tr.stats.npts >= tw.stats.npts:
  fct=xcorr(tr.data,tw.data)
  t,coef=xcor_cutTimes(tr,fct,xcor_cut,plot=True)

Reference:
Holland, A. A., 2013, Earthquakes Triggered by Hydraulic Fracturing in South-Central 
   Oklahoma: Bull. Seismol. Soc. Am., v. 103, no. 3, p. 1784-1792.

Author:
Austin Holland, austin.holland@ou.edu
Oklahoma Geological Survey, University of Oklahoma
2013,  
This work was funded by the Oklahoma Geological Survey.

This software is issued under the GPL http://www.gnu.org/licenses/gpl.txt.  
"""
from obspy.core import *
import numpy as np
from scipy.signal import argrelmax
import matplotlib.pyplot as plt

def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
    
def xcorr(x,y):
  """c=xcor(x,y)
  Fast implementation to compute the normalized cross correlation where x and y are 1D numpy arrays
  x is the timeseries
  y is the template time series
  returns a numpy 1D array of correlation coefficients, c"
  
  The standard deviation algorithm in numpy is the biggest slow down in this method.  
  The issue has been identified hopefully they make improvements.
  ""

  N=len(x)
  M=len(y)
  meany=np.mean(y)
  stdy=np.std(np.asarray(y))
  tmp=rolling_window(x,M)
  c=np.sum((y-meany)*(tmp-np.reshape(np.mean(tmp,-1),(N-M+1,1))),-1)/(M*np.std(tmp,-1)*stdy)

  return c

def xcor_cutTimes(trace,fct,cut,plot=False):
  """ xcor_cutTimes(trace,fct,cut,plot=False)
  Returns the times and maximum correlation values in an obspy trace (trace) above the 
  cut off value (cut) that maximize the correlation coefficients (fct) output from xcor.
  plot=True allows for visual confirmation of expected results
  """

#  t_indx=np.where(fct>=cut) # Where had some bad behaviors
  mfct=np.clip(fct,cut-.001,2)
  if np.max(mfct)==cut-.001:
    indx=[]
  else:
    max=argrelmax(mfct)
    indx=max[0]
  coef=[]


      #print cor
#  print indx,coef     
#  print cindx

#  print np.max(fct)
  if plot:
    plt.figure()
    plt.subplot(211)
    plt.plot(trace.data)
    plt.subplot(212)
    plt.plot(np.arange(0,len(fct)),fct)
    plt.plot(indx,np.ones(len(indx))*cut,'+r')
    plt.ylim((-1,1))
    plt.show()
  t=[]

  for ndx in indx:
    #print ndx
    t.append(trace.stats.starttime+(ndx*trace.stats.delta))
    coef.append(fct[ndx])
  return t,coef
  

