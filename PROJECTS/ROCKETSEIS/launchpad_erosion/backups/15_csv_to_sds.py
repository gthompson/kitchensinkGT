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

dryrun = False

sdsobj = SDS.SDSobj(paths['SDS_TOP'], sds_type='D', format='MSEED')

# Parse lookuptable and convert good CSV files to SDS
lookuptableDF = LLE.removed_unnamed_columns(pd.read_csv(paths['lookuptable']))
lookuptableDF.to_csv('lookuptable_backup.csv')
lookuptableDF = lookuptableDF.sort_values(by=['starttime'])
#print("lookuptableDF = \n", lookuptableDF.columns)
#print(lookuptableDF)


transducersDF = LLE.removed_unnamed_columns(pd.read_csv(paths['transducersCSVfile']))
#print("transducersDF = \n",transducersDF.columns)

lookuptableDF['SDS'] = False

sds_column_exists = False
if 'SDS' in lookuptableDF.columns:
    sds_column_exists = True

for index, row in lookuptableDF.iterrows():
    row_to_sds_done = False
    if row['passed']: # only convert to SDS the rows that passed the quality metrics
        successful = False
        if sds_column_exists:
            row_to_sds_done = row['SDS']
        print(f"{index}, {row['sourcefile']}, Passed? True, SDS done? {row_to_sds_done}")
        if not row_to_sds_done:
            if dryrun: # SCAFFOLD: do not process if dryrunning for now, because only interested in the files that did not pass
                continue 
            df2 = pd.read_csv(os.path.join(paths['CORRECTED'],row['outputfile']))
            successful = LLE.convert2sds(df2, sdsobj, transducersDF, dryrun=dryrun)
            lookuptableDF.at[index,'SDS'] = successful
            if successful:
                lookuptableDF.to_csv(paths['lookuptable']) 
    else: # SCAFFOLD: what to do with rows that failed to pass
        print(f"SCAFFOLD: {index}, {row['sourcefile']}, Passed? False, SDS done? False")      
        df2 = pd.read_csv(os.path.join(paths['CORRECTED'],row['outputfile']))
        successful = LLE.convert2sds_badfile(df2, sdsobj, transducersDF, dryrun=dryrun)    
        lookuptableDF.at[index,'SDS'] = successful
        if successful:
            lookuptableDF.to_csv(paths['lookuptable'])         
del sdsobj


