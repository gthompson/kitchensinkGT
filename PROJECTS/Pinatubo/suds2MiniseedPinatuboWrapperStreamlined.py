#!/usr/bin/env python
# coding: utf-8

# # SUDS to MiniSEED wrapper for all Pinatubo waveform data
# This Python3 notebook demonstrates how to convert a PC-SUDS seismic waveform data file into Miniseed.
# It used to require Win-SUDS Utilities (for programs demux.exe, irig.exe and sud2sac.exe), but ObsPy now 
# includes a SUDS DMX reader. After comparing this with data converted with sud2sac.exe, I wrapped this in my
# own function to match as closely as possible.

# Do all library imports first
import obspy.core as op
import glob
import matplotlib.pyplot as plt
import os
import shutil

# Set up all paths that will be used
os.chdir('D:\Dropbox\DATA\Pinatubo')
#os.putenv('SEISAN_TOP', 'D:\Dropbox\DATA\SEISAN_DB')
os.environ['SEISAN_TOP']='D:\Dropbox\DATA\SEISAN_DB'
net = 'XB' # assigned by Gale Cox. 1R is code for KSC.
seisanDBname = 'PINAT'
WORpath = os.path.join(os.getenv('SEISAN_TOP'),'WOR')
WAVEFORM_DIR = 'WAVEFORMS'    
print('Paths setup okay')


# Functions
def findLatestFile(dirpath):
    lists = os.listdir(dirpath)                                   
    lists.sort(key=lambda fn:os.path.getmtime(os.path.join(dirpath, fn)))
    return os.path.join(dirpath,lists[-1])

def displayFile(dirpath):
    if os.path.exists(dirpath):
        fptr = open(dirpath, 'r')
        str = fptr.read()
        print(str)
    else:
        print(dirpath + ' not found')
        
def fixNSLC(st):
    for tr in st:
        sta = tr.stats.station
        tr.stats.network = net
        tr.stats.station = sta[:-1]
        tr.stats.channel = 'EH%c' % sta[-1]
        tr.stats.sampling_rate=100.0
        print('Converting %s => %s' % (sta, tr.id))
        
def read_DMX_file(DMXfile):
    # I have found this produces same result as converting DMX to SAC with sud2sac.exe in Win-SUDS, and then reading into ObsPy
    # Tested only on Pinatubo 1991 data
    # DMXfile is normally read in as uint16 and so is all +ve. But SUDS data must actually be 2s-complement or sign bit, because
    # sud2sac.exe converts to numbers either size of 0. Difference always seems to be 2048. 
    # Obspy Miniseed writer needs float, not int, so will crash unless recast as float.
    # However, sampling rates between DMX and SAC are different
    st = op.read(DMXfile)
    for tr in st:
        if tr.stats.station=='IRIG':
            st.remove(tr) # we do not want to keep the IRIG trace
        if all(tr.data==0):
            st.remove(tr)
        if tr.stats.network == 'unk':
            tr.stats.network = ''
        tr.data = tr.data.astype(float) - 2048.0     
    return st

# Loop over all files
badDMXfiles = []
alldirs = glob.glob(os.path.join(WAVEFORM_DIR, '*'))
for thisdir in alldirs:
    allDMXfiles = glob.glob(os.path.join(thisdir, '*.DMX'))
    for DMXfile in allDMXfiles:
        print(DMXfile)
        DMXbasename = os.path.basename(DMXfile)
        try:
            st = op.read(DMXfile, headonly=True)
        except:
            print("- FAILED to read")
            badDMXfiles.append(DMXbasename)
            continue
        WAVbasename = "%sM.%s_%03d" % (st[0].stats.starttime.strftime('%Y-%m-%d-%H%M-%S'), seisanDBname, len(st))
        #PHAfile = SUDSbasename + '.PHA' # this might exist if HYPO71 was run to locate the event
        #PUNfile = SUDSbasename + '.PUN' # this might exist if HYPO71 was run and generated a hypocenter
        WAVfile = os.path.join(os.getenv('SEISAN_TOP'), 'WAV', seisanDBname, 
                               st[0].stats.starttime.strftime('%Y'), 
                               st[0].stats.starttime.strftime('%m'),
                               WAVbasename)
        if os.path.exists(WAVfile):
            print('%s exists' % WAVfile)
            continue

        # Read DMX file
        try:
            st = read_DMX_file(DMXfile)
            print('- read okay')
            print(st)
        except:
            print('- ObsPy cannot read this demultiplexed SUDS format file')
            continue
            
        if len(st)==0:
            print('Empty Stream object')
            badDMXfiles.append(DMXbasename)
            continue   
            
        # fix/map NSLC    
        fixNSLC(st)

        # Add this Miniseed event waveform file to the Seisan database.
        # The Miniseed file will be moved to [SEISANTOP]/WAV/[SEISANDB]/YYYY/MM
        # A corresponding Seisan S-file (event metadata file) will be created at [SEISANTOP]/REA/[SEISANDB]/YYYY/MM
        # The Seisan programs MAKEREA and AUTOREG are used here. Since they normally require user input, 
        # we create files to simulate this.  
        yy = DMXbasename[0:2]
        mm = DMXbasename[2:4]
        century = '19'
        if yy < '80':
            century = '20'
        yyyy = century + yy
        fptr = open('makerea_wrapper.txt','w')
        fptr.write(seisanDBname + '\n')
        fptr.write(yyyy + mm + '\n')
        fptr.write('\n')
        fptr.write('BOTH\n')
        fptr.close()
        os.system('makerea < makerea_wrapper.txt')

        # Copy the miniseed file to WOR BEFORE running dirf and autoreg
        WORfile = os.path.join(WORpath, WAVbasename)
        try: 
            st.write(WORfile,'mseed')
        except:
            print('FAILED to write')
            continue
        cwd=os.getcwd()
        print('cwd = ' + cwd)
        os.chdir(WORpath)
        os.system('dirf ' + os.path.basename(WORfile))
        fptr = open('autoreg_wrapper.txt','w')
        fptr.write('L\n')
        fptr.write('m\n')
        fptr.write(seisanDBname + '\n')
        fptr.write('gt\n')
        fptr.close()
        os.system('autoreg < autoreg_wrapper.txt')
        os.chdir(cwd)


      


# Can now browse event with:
# eev yyyymm PNTBO
# Can also view the latest S-File with the code below
SFILE = findLatestFile(os.path.join(os.environ['SEISAN_TOP'],'REA',seisanDBname,yyyy,mm))
print(SFILE)
displayFile(SFILE)  


# There are no *.PHA files with the Pinatubo dataset John Power has given me. Instead, there is:
# - pinmay91.pha, ..., pinaug91: these are probably concatenated PHA files from each month, but missing Sep-Dec.
# - PIN_STA_ALL.txt: station coordinates
# - Pinatubo_all.SUM: Summary file of all Pinatubo hypocenters
# - b-run.SUM: Perhaps an attempt to locate all events, some of which failed. Although a SUM file, different format to Pinatubo_all.SUM, and times of located events do not appear to correspond.
# 
# What remains to be done:
# 1. Some DMX files seem to be broken. Why?  process badDMXfileList.txt?
# 2. make the list of number of files each day from the last event label of each day
# 3. Add DOI's for dataset and software & create github project
# 4. Go through Pinatubo_all.SUM list, and try to associate with any event it overlaps with. Easiest to do this in ObsPy. 
# 5. Streamlining: Modify script so that we do not run makerea each time
# 6. Are DMX files stamped with local time, or UTC? 
# 7. Process pinmay91.pha, ..., pinaug91, try to match with existing S-Files.
# 8. Add station coordinates from  PIN_STA_ALL.txt to STATION0.HYP
# 9. Get instrument responses and add to Seisan.
