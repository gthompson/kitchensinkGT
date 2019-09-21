# Amplitude Source Location
# HISTORY:
# 20190921: initial template Glenn Thompson

#from obspy import read
##from obspy.core import Stream
#import os.path
#import numpy as np
#import datetime
#import sys
##import matplotlib.pyplot as plt
#import struct
#import os
#import glob
##import time
##import re
#import warnings
#from pathlib import Path

#warnings.filterwarnings("ignore")
##sys.path.insert(0,'~/scripts')
##import libMVOarchive

def get_stations():
    # stations and their:
    #     lat, lon, elev coordinates
    #     instrument corrections
    #     site corrections
    #     combined corrections (e.g. from regional earthquake amplitude analysis)
    return # stations dict

def get_next_WAV():
    # load an event from Seisan WAV file into a Stream object
    return # Stream object

def load_Sfile():
    # this should use methods from SEISAN/Seisan_Catalog.py
    return # Seisan_Catalog object

def initial_location():
    # choose an initial location, e.g. dome coordinates for Montserrat
    # perhaps vary depth based on event classification if available
    load_Sfile()
    return # location object

def update_location(previousloc, previouserrormetric, amp_observed, iteration_num):
    # normalize real amplitude distribution, amp_observed
    create_grid(previousloc, iteration_num) # create grid around loc with node spacing determined by zoomlevel
    # grid should shrink by irrational factor as iteration_num increases
    # set bestloc to previousloc
    # set minerrormetric to previouserrormetric
    # loop over each grid node
    #    set nodeloc to this node location
    #    compute expected normalized amplitude distribution at this node, amp_expected
    #    errormetric = distance between amp_observed and amp_expected
    #    if errormetric < minerrormetric, set minerrormetric to errormetric & bestloc to nodeloc
    return bestloc, minerrormetric # location object


def timewindow_locate(tw):
    # get initial location
    loc = initial_location()
    # reduce Stream object to mean/max/median/rms amplitude
    amp = stream2amplitude()
    num_iterations = 0
    errormetric = Inf
    errorthreshold = 1 # no idea what this should be
    while ((errormetric > errorthreshold) and (num_iterations < 10)):
        loc,errormetric = update_location(loc, errormetric, amp, num_iterations)
        num_iterations += 1
    return loc, errormetric
    
def get_next_timewindow(st, timewindow_start, timewindow_seconds):
    # extract next timewindow from Stream object 
    return # Stream object

def main(argv):
    wavfilepath = argv[1]
    get_stations()
    st = loadWAVfile(wavfilepath) 
    #stream_endtime = 
    # need to sort out different time formats here, probably not in seconds
    timewindow_seconds = 10
    timewindow_overlap = 5
    timewindow_start = 0
    while timewindow_start < stream_endtime:
        tw = get_next_timewindow(st, timewindow_start, timewindow_seconds)
        loc, errormetric = timewindow_locate(tw) # need to add these to array/list to track trajectory 
        timewindow_start = timewindow_start + timewindow_overlap
    

if __name__ == "__main__":
    # Usage PROGRAMNAME sourcedir seisan_top seisan_dbname
    main(sys.argv)
