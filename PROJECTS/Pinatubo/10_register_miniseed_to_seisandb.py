#!/usr/bin/env python
# coding: utf-8

import obspy.core as op
import glob
import os
#import shutil
#import numpy as np
import subprocess
#import pandas as pd
#import numpy as np
#from IPython.display import display, HTML

def _run_command(cmd, show_terminal_output=True, errorStr=''):
    success = True
    print(cmd)
    result = subprocess.run(cmd, stdout=subprocess.PIPE) #, check=True)  
    #print(result)
    if show_terminal_output:
        if result.stdout:
            #print('\n')
            print("STDOUT:", result.stdout.decode())  # decode the byte-string
            #print('\n')
        if result.stderr:
            print("STDERR:", result.stderr.decode())
            #print('\n')
    if result.stderr:
        return 1, result
    elif errorStr and result.stdout.decode().find(errorStr):
        return 2, result
    else:
        return 0, result
   
        
def fix_sampling_rate(st, fs=100.0):
    for tr in st:
        tr.stats.sampling_rate=fs          

def _Stream_to_SEISANWAV_filename_translation(st, SEISAN_TOP, seisanDBname):
    # Establishing the Seisan WAV filename corresponding to this Stream
    WAVbasename = "%sM.%s_%03d" % (st[0].stats.starttime.strftime('%Y-%m-%d-%H%M-%S'), seisanDBname, len(st))
    WAVfilename = os.path.join(SEISAN_TOP, 'WAV', seisanDBname, 
        st[0].stats.starttime.strftime('%Y'), 
        st[0].stats.starttime.strftime('%m'),
        WAVbasename)    
    return WAVfilename

def register_Stream_into_SeisanDB(st, SEISAN_TOP, seisanDBname):
    ### REGISTERING THIS EVENT FILE INTO SEISAN
    # Check if it already exists in Seisan DB WAV
    WAVfilename = _Stream_to_SEISANWAV_filename_translation(st, SEISAN_TOP, seisanDBname)
    
    if os.path.exists(WAVfilename):
        print('%s exists' % WAVfilename)
        return 1
    else:
        # Add this Miniseed event waveform file to the Seisan database.
        # The Miniseed file will be moved to [SEISANTOP]/WAV/[SEISANDB]/YYYY/MM
        # A corresponding Seisan S-file (event metadata file) will be created at [SEISANTOP]/REA/[SEISANDB]/YYYY/MM
        # The Seisan programs MAKEREA and AUTOREG are used here. Since they normally require user input, 
        # we create files to simulate this. 
        WAVbasename = os.path.basename(WAVfilename)
        #
        # MAKEREA
        # old stuff. was getting YYYY & MM from DMX file
        #DMXbasename = os.path.basename(originalDMXfile)
        #yy = DMXbasename[0:2]
        #mm = DMXbasename[2:4]
        #century = '19'
        #if yy < '80':
        #    century = '20'
        #yyyy = century + yy
        # New stuff. get YYYY & MM from starttime 
        yyyy = st[0].stats.starttime.strftime('%Y')
        mm = st[0].stats.starttime.strftime('%m')
        fptr = open('makerea_wrapper.txt','w')
        fptr.write(seisanDBname + '\n')
        fptr.write(yyyy + mm + '\n')
        fptr.write('\n')
        fptr.write('BOTH\n')
        fptr.close()
        cmd = 'makerea < makerea_wrapper.txt'
        rtncode, result = _run_command(cmd, show_terminal_output=True)
        
        #os.system('makerea < makerea_wrapper.txt')
        #
        # Copy the miniseed file to WOR BEFORE running dirf and autoreg
        WORpath = os.path.join(SEISAN_TOP,'WOR')
        WORfile = os.path.join(WORpath, WAVbasename)
        if not os.path.isfile(WORfile):
            try: 
                st.write(WORfile,'mseed')
            except:
                print('FAILED to write')
                return 2 
        #
        # DIRF    
        cwd=os.getcwd()
        print('cwd = ' + cwd)
        os.chdir(WORpath)
        #os.system('dirf ' + os.path.basename(WORfile))
        cmd = 'dirf ' + WAVbasename
        try:
            rtncode, result = _run_command(cmd, show_terminal_output=True)
            if rtncode:
                return 3
        except:
            return 4
        #
        # AUTOREG
        fptr = open('autoreg_wrapper.txt','w')
        fptr.write('L\n')
        fptr.write('m\n')
        fptr.write(seisanDBname + '\n')
        fptr.write('gt\n')
        fptr.close()
        #os.system('autoreg < autoreg_wrapper.txt')
        cmd = 'autoreg < autoreg_wrapper.txt'
        try:
            rtncode, result = _run_command(cmd, show_terminal_output=True)
            if rtncode:
                return 5
        except:
            return 6
        os.chdir(cwd)
        
        return 0


# Main paths
TOPDIR = 'D:\Dropbox\DATA\Pinatubo'
CONVERT_DIR = 'D:\convert'

# paths for files to register into Seisan
os.environ['SEISAN_TOP']='D:\Dropbox\DATA\SEISAN_DB'
SEISAN_TOP = os.getenv('SEISAN_TOP')
seisanDBname = 'PINAT'
#WORpath = os.path.join(os.getenv('SEISAN_TOP'),'WOR')


# Loop over all files
allWORfiles = glob.glob(os.path.join(CONVERT_DIR, '*M.%s*' % seisanDBname))
for worfile in allWORfiles:
    print('Processing %s' % worfile)
    st = op.read(worfile) 
    
    # Registering into Seisan
    rtncode = register_Stream_into_SeisanDB(st, SEISAN_TOP, seisanDBname)


# It does not seem that either makerea or autoreg work. So will have to run them from terminal. 
# In which case, this program is just a copy frm D:\convert to SEISAN_TOP/WOR
