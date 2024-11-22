#!/Users/thompsong/miniconda3/bin/python

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 13:02:50 2015

@author: thompsong
"""

from obspy import read
import os.path
import numpy as np
import datetime
import sys
import matplotlib.pyplot as plt
import struct
import os
import glob
#import scipy.stats as scistats

def summarize(arrays):
    eventdate, eventtime, duration, signalamp, noiseamp, snr = arrays
    total_duration = np.nansum(duration) 
    median_snr = np.nanmedian(snr)
    return (total_duration, median_snr, ) 

def get_filelist(yearmonthdir):
    files = sorted(glob.glob(os.path.join(yearmonthdir, 'MV.M*.*')))
    return files

def build_arrays(thisfile):
    f = open(thisfile)
    duration = np.array([])
    snr = np.array([])
    signalamp = np.array([])
    noiseamp = np.array([])
    for line in f:
        fields = line.strip().split()
        # Array indices start at 0 unlike AWK
        eventdate = fields[0]
        eventtime = fields[1]
        duration = np.append(duration, float(fields[2]))
        signalamp = np.append(signalamp, float(fields[3]))
        noiseamp = np.append(noiseamp, float(fields[4]))
        snr = np.append(snr, float(fields[5]))
    f.close()
    return (eventdate, eventtime, duration, signalamp, noiseamp, snr, )

def plot_signalamp(thisfile, signalamp):
    plt.plot(signalamp, '*')
    plt.ylabel('Signal amplitude (counts)')
    plt.xlabel('Event number')
    plt.title(thisfile)
    plt.show()
    
def plot_snr_vs_signalamp(thisfile, signalamp, snr):
    plt.plot(signalamp, snr, "*")
    plt.xlabel('Signal amplitude (counts)')
    plt.ylabel('Signal to noise ratio')
    plt.title(thisfile)
    plt.show()

def main():
    PROJECT_TOP = '/Users/thompsong/Desktop/IPGP_Thompson_collaboration'
    antelopetop = os.path.join(PROJECT_TOP, 'miniseed')
    seisantop = os.path.join(PROJECT_TOP, 'seismo')
    db = 'MVOE_' #seisan database code
    yearmonthdir = os.path.join(antelopetop, 'stats', '2005', '01')
    files = get_filelist(yearmonthdir)
    for thisfile in files:
        fname = os.path.basename(thisfile)
        arrays = build_arrays(thisfile)
        total_duration, median_snr = summarize(arrays)
        print('%s: \t%9.2f s \t%6.2f' % (thisfile, total_duration, median_snr, ) )
        eventdate, eventtime, duration, signalamp, noiseamp, snr = arrays
        #plot_signalamp(fname, signalamp)
        #plot_snr_vs_signalamp(fname, signalamp, snr)

if __name__ == "__main__":
	main()
