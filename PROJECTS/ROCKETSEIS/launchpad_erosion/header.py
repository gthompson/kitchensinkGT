#!/usr/bin/env python
# coding: utf-8

import os, sys

def setup_environment():
    paths = dict()
    paths['HOME'] = os.path.expanduser('~')
    paths['Developer'] = os.path.join(paths['HOME'], 'Developer')
    paths['tremor_explorer_lib'] = os.path.join(paths['Developer'], 'tremorExplorer', 'lib')
    paths['src'] = os.path.join(paths['Developer'], 'kitchensinkGT', 'PROJECTS', 'ROCKETSEIS', 'launchpad_erosion')
    sys.path.append(paths['tremor_explorer_lib'])
    sys.path.append(paths['src'])
    paths['work'] = os.path.join(paths['HOME'], 'work')
    paths['outdir'] = os.path.join(paths['work'], 'PROJECTS', 'KSC_EROSION')

    '''
    # ONLY NEED THIS PART FOR INITIAL CONVERSIONS FROM TOB3 FILES TO CSV (or PKL)
    # I FIRST USED CS python libraries, and then later CS LoggerNet. 
    # Do not need this anymore, and could just move raw TOB3 data to newton (and out of Dropbox)
    import platform
    if platform.system() == 'Windows':
        paths['DROPBOX_TOP'] = 'D:/Dropbox' 
    else:
        paths['DROPBOX_TOP'] = os.path.join(paths['HOME'], 'Dropbox')
    if os.path.isdir(paths['DROPBOX_TOP']):
        # For Well data conversion from Campbell Scientific datalogger format to Pickle Files, and eventually SDS
        paths['WELLDATA_TOP'] = os.path.join(paths['DROPBOX_TOP'], 'DATA', 'KSC', 'KSC_Well_Seismoacoustic_Data', 'WellData')
        paths['TOB3_DIR'] = os.path.join(paths['WELLDATA_TOP'], 'Uploads')
        sys.path.append(os.path.join(paths['src'], 'campbell'))
        #import read_cs_files as campbell
    '''
    
    # CORRECTED DATA - always need this
    paths['CORRECTED'] = os.path.join(paths['outdir'], 'corrected')
    paths['lookuptable'] = os.path.join(paths['CORRECTED'],'lookuptable.csv') 
    if not os.path.isdir(paths['CORRECTED']):
        os.mkdir(paths['CORRECTED'])
    sys.path.append(os.path.join(paths['src'], 'lib'))
    # CORRECTED DATA is in PSI. Need to convert to depth in feet when plotting.
    paths['transducersCSVfile'] = os.path.join(paths['outdir'], 'transducer_metadata.csv')

    '''
    # read general iceweb config
    paths['CONFIGDIR'] = os.path.join(paths['Developer'], 'tremorExplorer', 'config')
    if not 'PRODUCTS_TOP' in locals():
        PRODUCTS_TOP=None
    configDict = wrappers.read_config(configdir=paths['CONFIGDIR'], PRODUCTS_TOP=PRODUCTS_TOP)
    paths['SDS_TOP'] = configDict['general']['SDS_TOP']
    paths['RSAM_DIR'] = configDict['general']['RSAM_TOP']
    '''
    paths['SDS_TOP'] = os.path.join(paths['outdir'], 'SDS')
    paths['DB_DIR'] = os.path.join(paths['outdir'], 'db')
    paths['DB_DIR'] = os.path.join(paths['outdir'], 'RSAM')
    return paths

if __name__ == "__main__":
    paths = setup_environment()
    print(paths)
    print("sys.path\n",sys.path)
    print('Testing imports')
    import FDSNtools
    import wrappers
    import SDS
    import libWellData as LLE    

    
# These are additional imports to make
#import os
#import sys
#import glob
#import numpy as np
#import pandas as pd
#from obspy.core import read, Stream, UTCDateTime
#import FDSNtools
#import wrappers
#import SDS
#import libWellData as LLE

