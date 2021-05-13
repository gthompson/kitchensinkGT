#!/usr/bin/env python
# coding: utf-8

# # IRIS PASSCAL - RT130 Data Processing
# 
# Note: this version was built using a 2020 version of rt2ms, which is significantly streamlined compared to the workflow with 2016 version of rt2ms. The two main advantages are (1) rt2ms can build the parameter file, rather than constructing a dbbuild_batch.pf file and using batch2par to convert it, and (2) rt2ms directly builds daily MSEED volumes in BUD format, rather than building hourly MSEED files which then had to be combined into daily MSEED using dataselect.
# 
# You’ve offloaded a service run and have data from each RT130. Use this program to convert the data to miniSEED and reorganize it into station/channel/day volumes. Then, create a stationXML for your experiment using Nexus before submitting data to PASSCAL. Additional documentation can be found on the PASSCAL website: https://www.passcal.nmt.edu/content/passive-source-seed-archiving-documentation

# Import needed modules and set global variables

import os
#import glob
#import shutil
#from pathlib import Path
#import numpy as np
#import pandas as pd
print('Module imports - Done') 

network = '1R'

# ## STEP 1. Create an organized directory structure for your data. Below is a snippet of the data organization assumed. The variable $REFTEKDIR should point to the ROCKETSEIS directory in this example. 
"""
ROCKETSEIS/
└── SVC2
    └── RAW
        ├── 2018289
        │   ├── 91F8
        │   │   ├── 0
        │   │   ├── 1
        │   │   └── 9
        │   ├── 9338
        │   │   ├── 0
        │   │   ├── 1
        │   │   └── 9
        │   ├── 9406
        │   │   ├── 0
        │   │   ├── 1
        │   │   └── 9
        │   ├── 98E6
        │   │   ├── 0
        │   │   ├── 1
        │   │   └── 9
        │   ├── 9BC5
        │   │   ├── 0
        │   │   ├── 1
        │   │   └── 9
        │   ├── 9D7C
        │   │   ├── 0
        │   │   ├── 1
        │   │   └── 9
        │   └── AB13
        │       ├── 0
        │       ├── 1
        │       └── 9
        └── 2018290
            ├── 91F8
            │   ├── 0
            │   ├── 1
            │   └── 9
            ├── 9338
            │   ├── 0
            │   ├── 1
            │   └── 9
            ├── 9406
            │   ├── 0
            │   ├── 1
            │   └── 9
            ├── 98E6
            │   ├── 0
            │   ├── 1
            │   └── 9
            ├── 9BC5
            │   ├── 0
            │   ├── 1
            │   └── 9
            ├── 9D7C
            │   ├── 0
            │   ├── 1
            │   └── 9
            └── AB13
                ├── 0
                ├── 1
                └── 9
"""
# Here SVC2 stands for service run 2, and there is a directory RAW underneath this, and under that we have 7 4-character digitizer directories, and each of those has 0, 1 and 9 directories. 
# To create this directory structure, rather than using *Neo* to read the Compact Flash cards, which has the disadvantage of compressing the data to ZIP format (which is time consuming), and then saving the ZIP file to the computer, I just copied the data directory at the bash command line with:
# 
#     cp -r /Volumes/UNTITLED/RT130-*/2* ROCKETSEIS/SVC2/RAW/


REFTEKDIR = '/raid/data/KennedySpaceCenter/duringPASSCAL/REFTEK_DATA/ROCKETSEIS' 
SERVICE_RUN = 2 

yn = input('Is %s path to the main project directory correct? ' % REFTEKDIR)
if not (yn.lower() == 'yes' or yn.lower() == 'y'):
    REFTEKDIR = input('Enter correct path to main project directory: ')
    
yn = input('Is this service run %d ? ' % SERVICE_RUN)
if not (yn.lower() == 'yes' or yn.lower() == 'y'):
    SERVICE_RUN = int(input('Enter correct service run number: '))
       
print('Setting paths for relative directories/files')
SVCDIR = os.path.join(REFTEKDIR, 'SVC%d' % SERVICE_RUN)
RAWDIR = os.path.join(SVCDIR, 'RAW') 
#LOGSDIR =  os.path.join(SVCDIR, 'LOGS')
#CONFIGDIR = os.path.join(SVCDIR, 'CONFIG')
#MSEEDDIR = os.path.join(SVCDIR, 'MSEED')
print('Path to raw data is %s\n' % RAWDIR)

# *commandExists* checks if PASSOFT/DMC commands are installed before we try to use them.
def commandExists(command):
    output = os.popen('which %s' % command).read()
    if output:
        return True
    else:
        print('Command %s not found.' % command)
        print('Make sure the PASSOFT tools are installed on this computer, and available on the $PATH')
        return False


# ### STEP 3: Convert your data into miniSEED files. 
# In the service run directory, convert the raw RT130 data to miniSEED. Typing *rt2ms -h* shows a list of available options. 
yn = input('Are you ready to run rt2ms ? ')
if not (yn.lower() == 'yes' or yn.lower() == 'y'):
    raise SystemExit("Killed!") 
parfile = 'parfile.txt' # this name is fixed by rt2ms
if os.path.exists(parfile): 
    os.system('rm %s' % parfile) 
os.system('ls -d %s/20????? > dir.list' % RAWDIR)
os.system('rt2ms -D dir.list -e') # creates a file called parfile.txt
if os.path.exists(parfile): 
    # Edit the par file
    os.system("sed -i -e 's/XX/%s/' %s" % (network, parfile))
    os.system("sed -i -e 's/ELZ/EHZ/' %s" % parfile)
    os.system("sed -i -e 's/EL1/EH1/' %s" % parfile)
    os.system("sed -i -e 's/EL2/EH2/' %s" % parfile)
    os.system('rt2ms -D dir.list -p %s' % parfile)


# ## STEP 4: Confirm your station and channel names
# In the DAYS folder just created by dataselect, check to see if you have folders for each of your stations. The data should be organized into those folders in station/channel/day volumes named STA.NET.LOC.CHAN.YEAR.JULDAY. For example: BA01.XR..HHZ.2018.039 (The .. after XR is where the location code would be if needed).
# 
# If your parfile was incomplete (i.e. missing stations or channels), there will be one or more folders named with the RT130 serial numbers (e.g. 9306) instead of the desired station name (e.g. ME42). To change any miniSEED headers to correct a station name, network code, etc., see the *fixhdr* doc on the PASSCAL website (see link on the first page). After you have modified the headers with *fixhdr*, rename the files so that the station‐network‐location‐channel codes in the miniSEED file names match the corrected headers.  

print('Step 4: Confirm trace ids / BUD format. Use fixhdr to correct.\n')

# ## STEP 5. Perform quality control of waveforms and logs. 
# Verify the data quality by reviewing the traces and log files (with *logpeek* and *pql*). Obvious signs of trouble include loss of GPS timing, overlaps, gaps, corrupted files, etc. Make a note of any problems. Use *fixhdr* to correct mark timing issues, and/or to convert the files to big endianess if they are not already. For more information on how to use these tools, refer to the appropriate documentation on the PASSCAL website (see link on the first page).

print('Step 5: QC the waveform files. PASSCAL recommends using logpeek for log file QC and pql for waveform data. But I could not get logpeek to work. And easier towrap the dayfiles with an Antelope database and then use dbpick. Then use fixhdr to mark any problems.\n')

# ## STEP 6. Create metadata for your experiment. 
# Use *Nexus* to generate a stationXML file for your experiment metadata. See the “Metadata Generation with Nexus in a Nutshell” document on the PASSCAL website (see link on first page).
# 
# Essentially, you open Nexus & scan a set of Miniseed files. You have to manually enter the datalogger and sensor details, and the coordinates, depth etc. And any sensor misorientation information. Then compute responses. And save the stationXML file, e.g. to 1R.LOC00.xml.
print('Step 6: Use Nexus to create your stationXML file\n')

# ## STEP 7. Send miniSEED data to PASSCAL. 
# Please drop a note, with your PASSCAL project name in the subject, to <mailto>data_group@passcal.nmt.edu</mailto> before sending the data to PASSCAL so that we can set up a receiving area. Attach the stationXML created with Nexus to this email unless it is larger than 5Mb. Use our tool *data2passcal* to send the data: 
# 
#     data2passcal DAYS/ 
# 
# *data2passcal* will scan all subdirectories of the DAYS folder and send any miniSEED files that have the correct file names.
# 
print('Step 7: Use data2passcal to send DAYS directory to PASSCAL.\n')
