#!/usr/bin/env python
# coding: utf-8

import header
paths = header.setup_environment()
import os
#import sys
import glob
#import numpy as np
#import pandas as pd
#from obspy.core import read, Stream, UTCDateTime
#import FDSNtools
#import wrappers
#import SDS
#import libWellData as LLE

os.chdir(paths['outdir'])
import datetime

#ANTELOPE = os.getenv('ANTELOPE')
years = ['2022']
nets = ['6I', '6S', 'FL', 'AM', 'XA']
YYYY = '2022'

def dbcreate(dbname, schema="css3.0", dbpath=None):
    nl = '\n'
    dbpaths = ""
    if dbpath:
        for p in dbpath:
            d = os.path.dirname(p)
            b = os.path.basename(p)
            if dbpaths:
                connector = ":"
            else:
                connector = ""
            #thisdbpath = "%s/{%s}" % (d,b) # no, in same directory as descriptor, so just need dbbasename
            thisdbpath = "{"+b+"}"
            dbpaths += connector + thisdbpath
    contents = "#\n"
    contents += f"schema {schema}{nl}"
    if dbpath:
        contents += f"dbpath {dbpaths}{nl}"
    print(contents)
    with open(dbname, "w") as f:
        f.write(contents)

for YYYY in years:
    for jday in range(83,367): # 83 is first day of seismic data in SDS 
        start_date = datetime.datetime.strptime(f"{YYYY}{jday:03d}", '%Y%j')
        print(start_date, end="\n")
        ymd = start_date.strftime("%Y%m%d")
        dbday = f"db/db{ymd}"
        for net in nets:
            mseed_dirs = sorted(glob.glob(os.path.join('SDS', YYYY, net, '???*', '[HD]??.D')))
            for mseed_dir in mseed_dirs:
                mseedfiles = sorted(glob.glob(os.path.join(mseed_dir, '%s.*.D.*.%03d' % (net, jday) ) ))    
                if len(mseedfiles)>0:
                    mseedfilesstr = " ".join(mseedfiles)
                    os.system(f"miniseed2db {mseedfilesstr} {dbday} ")
        # create descriptor
        if os.path.isfile(dbday+'.wfdisc'):
            dbcreate(dbday, schema='css3.0')  


