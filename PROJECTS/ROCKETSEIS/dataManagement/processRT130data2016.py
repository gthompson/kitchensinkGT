#!/usr/bin/env python
# coding: utf-8

# # IRIS PASSCAL - RT130 Data Processing
# 
# A Jupyter notebook by Glenn Thompson based on: https://www.passcal.nmt.edu/webfm_send/3035
#
# Note: this version was built using a 2016 version of rt2ms. There is a more streamlined version for the 2020+ version of rt2ms.
# 
# You’ve offloaded a service run and have data from each RT130. Follow the steps in this document to convert the data to miniSEED and reorganize it into station/channel/day volumes. Then, create a stationXML for your experiment using Nexus (see step 7) before submitting data to PASSCAL. Program names are in italics. Unix commands and any command line arguments are on separate lines. Input files are denoted by < filename>. Additional documentation can be found on the PASSCAL website: https://www.passcal.nmt.edu/content/passive-source-seed-archiving-documentation

# Import needed modules and set global variables

import os
import glob
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
print('Done') 


# ## STEP 1. Create an organized directory structure for your data. 
# Start by creating a main directory for the project *(in this Jupyter Notebook, I use the variable REFTEKDIR for this main project directory)*. 
# 
# Under your main project directory, make a first level directory “SVC1” for service run number 1. For each subsequent service run create a new directory, e.g. SVC2, SVC3. Create directories in the SVC1 directory for the raw data files and log files. For example: 
# 
#     mkdir RAW
#     mkdir LOGS
# 
# Move the raw data files (either .ZIP or CF folders) into the RAW directory, e.g.
# 
#     mv SVC1/ZIPFILES/*.ZIP SVC1/RAW/
# 
# <b>Glenn's variations:</b> 
# <ol>
# <li>Rather than creating a directory for each service run, I created a new directory anytime there was a network change, e.g. station installed, sensor reoriented, DAS swapped etc. So these are directories like LOC00, LOC01, LOC02, etc. rather than SVC1, SVC2, SVC3.</li>
#  
# <li>While I sometimes used *Neo* to read the Compact Flash cards, compress the data to ZIP format, and then copy the data to the laptop, I mostly just copied the data using the MacOS command line as I found it quicker, e.g.</li>
# 
#     cp -r /Volumes/UNTITLED/RT130-*/2* LOC00/RAW/
#     
# <li>At the time of writing this Jupyter Notebook, I had already completed the field project and had all the data organized into a directory structure that looks like:
# </li></ol>
"""
LOC00/
    RAW/
        2018288/
            AB13/
                0/
                1/
                9/
            9D7C/
                0/
                1/
                9/
            ...
        2018289/
                ...
        ... 
LOC01/
    RAW/
        ...
"""
# <b>Note:</b> The REFTEKDIR variable below must be set to the path of the parent directory of SVC1 (or SVC2, SVC3, etc.). If you set SERVICE_RUN = 1, the Jupyter Notebook will create the SVC1 directory underneath this (or SVC2 if you set it to 2).
# 
# Note that I have modified this to map to LOC00, LOC01 etc., rather than SVC1, SVC2 etc.



REFTEKDIR = '/raid/data/KennedySpaceCenter/duringPASSCAL/REFTEK_DATA' 
SERVICE_RUN = 1

yn = input('Is %s path to the main project directory correct? ' % REFTEKDIR)
if not (yn.lower() == 'yes' or yn.lower() == 'y'):
    REFTEKDIR = input('Enter correct path to main project directory: ')
    
yn = input('Is this service run %d ? ' % SERVICE_RUN)
if not (yn.lower() == 'yes' or yn.lower() == 'y'):
    SERVICE_RUN = int(input('Enter correct service run number: '))
       
print('Setting paths for relative directories/files')
#SVCDIR = os.path.join(REFTEKDIR, 'SVC%d' % SERVICE_RUN)
SVCDIR = os.path.join(REFTEKDIR, 'LOC%02d' % (SERVICE_RUN - 1))
RAWDIR = os.path.join(SVCDIR, 'RAW') 
LOGSDIR =  os.path.join(SVCDIR, 'LOGS')
CONFIGDIR = os.path.join(SVCDIR, 'CONFIG')
MSEEDDIR = os.path.join(REFTEKDIR, 'MSEED')
DAYSDIR = os.path.join(REFTEKDIR, 'DAYS')
RT2MS_OUTPUT = os.path.join(LOGSDIR, 'rt2ms.out')

print('Creating outline directory structure')
if not os.path.exists(REFTEKDIR):
    print("%s does not exist. Exiting" % REFTEKDIR)
    raise SystemExit("Killed!")  
subdirs = [SVCDIR, RAWDIR, LOGSDIR, CONFIGDIR, MSEEDDIR, DAYSDIR]
for thissubdir in subdirs:
    if not os.path.exists(thissubdir):
        print('Need to make %s' % thissubdir)
        os.mkdir(thissubdir)
        if not os.path.exists(thissubdir):
            print("%s does not exist & could not be created. Exiting" % thissubdir)
            raise SystemExit("Killed!") 
    else:
        print('%s exists' % thissubdir)
yn = input('Do these paths look correct ? ')
if not (yn.lower() == 'yes' or yn.lower() == 'y'):
    raise SystemExit("Killed!") 

if os.path.exists(RT2MS_OUTPUT):
    os.remove(RT2MS_OUTPUT)
print('Outline directory structure created/exists')
print('Done')


# ## STEP 2. Create the parameter file(s). 
# The parameter file is used by *rt2ms* to assign header information to the miniSEED files. *rt2ms* is a PASSCAL program that gene
# rates miniSEED formatted files from REFTEK RT130 raw files. In addition, *rt2ms* also modifies the headers. In the SVC1 directory, use a text editor and information from your field notes to create an ASCII parameter file (parfile) following the examples at https://www.passcal.nmt.edu/webfm_send/3035.
# 
# <b>Glenn's variations</b>: PASSCAL instructions assume you construct a par file by hand for each network layout. However, I construct an Antelope-style *dbbuild_batch* pf file by hand for each network layout, and then use the PASSOFT program *batch2par* to convert this to a par file (in combination with 2 *sed* (Unix stream editor) commands to fix this. You can see *batch2par* and *sed* commands used below.


print('For each dbbuild_batch pf file, create a corresponding par file')

def chooseCorrectParFile(yyyyjjj):
    correctparfile = ""
    parfilelist = sorted(glob.glob('%s/network20?????.par' % CONFIGDIR))
    for thisparfile in parfilelist:
        parfileyyyyjjj = thisparfile[-11:-4]
        #print(parfileyyyyjjj, yyyyjjj)
        if parfileyyyyjjj <= yyyyjjj:
            correctparfile = thisparfile
    return correctparfile

def commandExists(command):
    output = os.popen('which %s' % command).read()
    if output:
        return True
    else:
        print('Command %s not found.' % command)
        print('Make sure the PASSOFT tools are installed on this computer, and available on the $PATH')
        return False

#pffiles = sorted(glob.glob('%s/network*.pf' % CONFIGDIR))
pffiles = sorted(glob.glob('%s/locations*.pf' % CONFIGDIR))
for pffile in pffiles:
    parfile = pffile[:-2] + 'par' # 
    print('- %s, %s' % (pffile, parfile)) 
            
    # Create the corresponding parfile if it does not already exist
    if not os.path.exists(parfile): 
        if commandExists('batch2par'): 

            print("batch2par %s -m > %s\n" % (pffile, parfile))
            os.system("batch2par %s -m > %s" % (pffile, parfile))            
            if os.path.exists(parfile): 
                # Edit the par file
                os.system("sed -i -e 's/rs200spsrs;/1;         /g' %s" % parfile)
                os.system("sed -i -e 's/x1/32/g' %s" % parfile);  
                os.system("sed -i -e '$s/None//' %s" % parfile);  
            else:
                print("- batch2par failed")
                raise SystemExit("Killed!")
        else:
            raise SystemExit("Killed!")
print('All pf files now have a corresponding par file')
print('Done')


# ### helper functions for part 3
# 
# *chooseCorrectParFile* selects the par file with YYYYJJJ closest to, but not exceeding, the day of Reftek data we are trying to process. For example, for 2018200 (200th day of year 2018), we match network2018200.par if it exists. Otherwise, network2018199.par is the best match, if it exists. Or whatever the closest pf is *before* the data day. But network2018201.par is not a match, since the existence of that pf would indicate the network changed on that day, so it was no longer applicable to 2018200 data.
# 
# *commandExists* checks if PASSOFT/DMC commands are installed before we try to use them.

# ### STEP 3: Convert your data into miniSEED files. 
# In the service run directory, convert the raw RT130 data to miniSEED. Typing *rt2ms -h* shows a list of available options. 
# 
# If raw data is in decompressed folders, use the following commands: 
# 
#     ls -d SVC1/RAW/*.cf > file.lst rt2ms -F file.lst -Y -L -o MSEED/ -p <parfile> >& rt2ms.out 
# 
# The (-F) flag will process all files in the named list, (‐Y) puts the data in yearly directories, (-L) outputs .log and, if created, .err files, (‐o) creates an output directory, MSEED, and (‐p) points to your parfile. 
# 
# If raw data is in ZIP files: 
# 
#     rt2ms ‐D SVC1/RAW/ ‐Y ‐L -o MSEED/ ‐p <parfile> >& rt2ms.out 
# 
# The (‐D) flag will process all .ZIP files in a specified directory, instead of in a file list as in the previous example. 
# 
# When *rt2ms* finishes, move all of your log and .err files from the MSEED directory to the LOGS directory that you created in step 1. 
# 
# After running *rt2ms* the MSEED directory structure should look something like the example below. In the MSEED directory there will be .log files and possibly .err files along with a  subdirectory for each year that contains day directories for each stream.
"""
MSEED/
    2014.019.21.29.16.98EZ.log
    2014.019.21.29.16.98EZ.err
    Y2014/
        R065.01/
        R065.02/
        R065.03/
"""
yn = input('Are you ready to run rt2ms ? ')
if not (yn.lower() == 'yes' or yn.lower() == 'y'):
    raise SystemExit("Killed!") 

def chooseCorrectParFile(yyyyjjj):
    correctparfile = ""
    #parfilelist = sorted(glob.glob('%s/network20?????.par' % CONFIGDIR))
    parfilelist = sorted(glob.glob('%s/locations20*.par' % CONFIGDIR))
    for thisparfile in parfilelist:
        parfileyyyyjjj = thisparfile[-11:-4]
        #print(parfileyyyyjjj, yyyyjjj)
        if parfileyyyyjjj <= yyyyjjj:
            correctparfile = thisparfile
    return correctparfile


dayfullpaths = sorted(glob.glob('%s/20?????' % RAWDIR))
for thisdayfullpath in dayfullpaths:
    print('Processing %s' % thisdayfullpath)
    thisdaydir = os.path.basename(thisdayfullpath) # a directory like 2018365
    
    # Find the corresponding parameter file for this day of the experiment
    parfile = chooseCorrectParFile(thisdaydir)
    if os.path.exists(parfile):                
        # Par file must exist if we got here. Run rt2ms.
        # Here we would ideally check that we have a full set of corresponding 0/, 1/ and 9/ files,
        # and we would also check if the hourly MSEED file already exists, and only run this if it does not
        if commandExists('rt2ms'): 
            rt130files = sorted(glob.glob('%s/????/[19]/*' % thisdayfullpath))
            if rt130files:

                for rt130file in rt130files:
                    try:
                        # We can check the following RT2MS_OUTFILE if rt2ms fails
                        # process one file at a time
                    	os.system("rt2ms -f %s -Y -L -p %s -o %s >> %s" % (rt130file, parfile, MSEEDDIR, RT2MS_OUTPUT))
                    except:
                        os.system("echo %s >> %s/badreftekfiles.txt" % (rt130file, REFTEKDIR))	

                # move all *.log files to the LOGS directory
                for src_file in Path(MSEEDDIR).glob('*.log'):
                    shutil.copy(src_file, LOGSDIR)

                # move all *.err files to the LOGS directory
                for src_file in Path(MSEEDDIR).glob('*.err'):
                    shutil.copy(src_file, LOGSDIR)         
                
        else:
            raise SystemExit("Killed!")
        
    else:
        print('- no corresponding parameter file found')
        
if dayfullpaths:
    print('FINISHED CONVERTING REFTEK DATA TO MINISEED HOURLY FILES')
else:
    print('No directories like RAW/YYYYJJJ found')
print('Done')    


# ## STEP 4: Reorganize the miniSEED data into station/channel/day volumes.
# 
# *dataselect* is a DMC program that allows for the extracting and sorting of miniSEED data (https://github.com/iris-edu/dataselect). This will read the data from the MSEED directory and convert them into day volumes with the required naming format: 
# 
#     dataselect -A DAYS/%s/%s.%n.%l.%c.%Y.%j MSEED/Y*/*/* 
# 
# The (-A) flag writes file names in the specified custom format. The format flags are (s) for station, (n) for netcode, (l) for location, (c) for channel name, (Y) for year, and (j) for Julian date. See the help menu for more details on options (*dataselect -h*). Depending on how much data you have, you may need to run *dataselect* in a loop that runs over the different days or stations in your experiment.
# 
# <b>Please note:</b> PASSCAL want data to be organized in BUD format. *dataselect -h* reveals that there is a BUD format option built in directly, so I attempt to use that here instead. This command is easier: 
# 
#     dataselect -BUD MSEED/Y*/*/* 
#     
# However, I do loop over year and day directories (as suggested) so that *dataselect* is not trying to deal with a file list that is too long for it or the operating system to handle.


yn = input('Are you ready to run dataselect to merge the hourly miniseed files into day miniseed files ? ')
if not (yn.lower() == 'yes' or yn.lower() == 'y'):
    raise SystemExit("Killed!") 

"""if commandExists('dataselect'): 
    # SCAFFOLD. I need to install this from github. And then figure out how to
    # substitute for station, netcode, location, channel name, year and julian day
    # May need to lop over different stations (I already do by day)
    mseedyearfullpaths = sorted(glob.glob('%s/Y20??' % MSEEDDIR))
    for thismseedyearfullpath in mseedyearfullpaths:
        # We only want to process dates that correspond to dayfullpaths

        print('Processing %s' % thismseedyearfullpath)  
        mseeddayfullpaths = sorted(glob.glob('%s/R*.01' % thismseedyearfullpath))
        for thismseeddayfullpath in mseeddayfullpaths:
            print('dataselect: Processing %s' % thismseeddayfullpath)
            #os.system('dataselect  -A %s/%s/%s.%n.%l.%c.%Y.%j  %s/*.m' % (DAYSDIR, thismseeddayfullpath))
            os.system('dataselect -BUD %s %s/*.m' % (DAYSDIR, thismseeddayfullpath) )
print('Done')  
"""           
if commandExists('dataselect'): 
    dayfullpaths = sorted(glob.glob('%s/20?????' % RAWDIR))
    for thisdayfullpath in dayfullpaths:
        yyyyjjj = os.path.basename(thisdayfullpath)
        thismseeddayfullpath = '%s/Y%s/R%s.01' % (MSEEDDIR, yyyyjjj[0:4], yyyyjjj[4:],  ) 
        print('dataselect: Processing %s' % thismseeddayfullpath)
        #os.system('dataselect  -A %s/%s/%s.%n.%l.%c.%Y.%j  %s/*.m' % (DAYSDIR, thismseeddayfullpath))
        os.system('dataselect -v -BUD %s %s/*.m' % (DAYSDIR, thismseeddayfullpath) )
print('Done')             


# ## STEP 5: Confirm your station and channel names
# In the DAYS folder just created by dataselect, check to see if you have folders for each of your stations. The data should be organized into those folders in station/channel/day volumes named STA.NET.LOC.CHAN.YEAR.JULDAY. For example: BA01.XR..HHZ.2018.039 (The .. after XR is where the location code would be if needed).
# 
# If your parfile was incomplete (i.e. missing stations or channels), there will be one or more folders named with the RT130 serial numbers (e.g. 9306) instead of the desired station name (e.g. ME42). To change any miniSEED headers to correct a station name, network code, etc., see the *fixhdr* doc on the PASSCAL website (see link on the first page). After you have modified the headers with *fixhdr*, rename the files so that the station‐network‐location‐channel codes in the miniSEED file names match the corrected headers.  

print('Step 5: Confirm trace ids / BUD format. Use fixhdr to correct.\n')
# ## STEP 6. Perform quality control of waveforms and logs. 
# Verify the data quality by reviewing the traces and log files (with *logpeek* and *pql*). Obvious signs of trouble include loss of GPS timing, overlaps, gaps, corrupted files, etc. Make a note of any problems. Use *fixhdr* to correct mark timing issues, and/or to convert the files to big endianess if they are not already. For more information on how to use these tools, refer to the appropriate documentation on the PASSCAL website (see link on the first page).

print('Step 6: QC the waveform files. PASSCAL recommends using logpeek for log file QC and pql for waveform data. But I could not get logpeek to work. And easier towrap the dayfiles with an Antelope database and then use dbpick. Then use fixhdr to mark any problems.\n')

# ## STEP 7. Create metadata for your experiment. 
# Use *Nexus* to generate a stationXML file for your experiment metadata. See the “Metadata Generation with Nexus in a Nutshell” document on the PASSCAL website (see link on first page).
# 
# Essentially, you open Nexus & scan a set of Miniseed files. You have to manually enter the datalogger and sensor details, and the coordinates, depth etc. And any sensor misorientation information. Then compute responses. And save the stationXML file, e.g. to 1R.LOC00.xml.
print('Step 7: Use Nexus to create your stationXML file\n')

# ## STEP 8. Send miniSEED data to PASSCAL. 
# Please drop a note, with your PASSCAL project name in the subject, to <mailto>data_group@passcal.nmt.edu</mailto> before sending the data to PASSCAL so that we can set up a receiving area. Attach the stationXML created with Nexus to this email unless it is larger than 5Mb. Use our tool *data2passcal* to send the data: 
# 
#     data2passcal DAYS/ 
# 
# *data2passcal* will scan all subdirectories of the DAYS folder and send any miniSEED files that have the correct file names.
# 

print('Step 8: Use data2passcal to send DAYS directory to PASSCAL.\n')




