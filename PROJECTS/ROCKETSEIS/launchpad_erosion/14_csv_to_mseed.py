#!/usr/bin/env python
# coding: utf-8
import header
paths = header.setup_environment()
import os
#import sys
#import glob
#import numpy as np
import pandas as pd
#from obspy.core import read, Stream, UTCDateTime
#import FDSNtools
#import wrappers
import SDS
import libWellData as LLE

# Parse lookuptable
lookuptableDF = LLE.removed_unnamed_columns(pd.read_csv(paths['lookuptable']))
lookuptableDF.to_csv('lookuptable_backup.csv')
lookuptableDF = lookuptableDF.sort_values(by=['starttime'])
lookuptableDF['miniseed'] = False

transducersDF = LLE.removed_unnamed_columns(pd.read_csv(paths['transducersCSVfile']))
MSEED_DIR = os.path.join(paths['outdir'], 'miniseed')
os.system(f"rm -rf {MSEED_DIR}/*")
for index, row in lookuptableDF.iterrows():
    print(f"{index}, {row['sourcefile']}, {row['passed']}")      
    df2 = pd.read_csv(os.path.join(paths['CORRECTED'],row['outputfile']))
    if row['passed']: # only convert to SDS the rows that passed the quality metrics
        successful = LLE.convert2mseed(df2, MSEED_DIR, transducersDF)
    else: # SCAFFOLD: what to do with rows that failed to pass
        successful = LLE.convert2mseed_badfile(df2, MSEED_DIR, transducersDF)    
    lookuptableDF.at[index, 'miniseed'] = successful
lookuptableDF.to_csv(paths['lookuptable'])   

