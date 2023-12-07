#!/usr/bin/env python
'''
Wrapper to walk Seisan WAV/DB/YYYY/MM tree and call autoreg on each WAV file
Glenn Thompson 2021/12/24
Caveats:
    1. Does not check if WAV file already registered
    2. Does not check if WAV file is compressed with gz or bzip2
Both checks could be easily added
'''
import os, glob

def autoreg_this_dir(thisdir, yyyy, autoreg_infile, fbad, n_good):
    os.chdir(thisdir) # otherwise S-file records absolute path
    # Loop over each WAV file
    wavfiles = glob.glob(os.path.join(thisdir, '%s*S.*' % yyyy )) # THIS FILE PATTERN MUST MATCH YOUR WAV FILES
    print('Processing %s: %d WAV files found' % (thisdir, len(wavfiles)))
    for wavfile in sorted(wavfiles):
        wavbase = os.path.basename(wavfile) 
        os.system('dirf %s' % wavbase) # Run dirf - produces filenr.lis with only 1 event
        try: # run autoreg
            os.system('autoreg < %s' % autoreg_infile)
            print('Registered %s' % wavbase)
            n_good += 1
        except: # it crashed. add to list of bad WAV files
            print('autoreg crashed on %s' % wavbase)
            fbad.write('%s\n' % wavfile)
    print('registered %d WAV files so far' % n_good)


# main program starts here

#--- Change these variables as needed ---#
seisantopdir = os.getenv('SEISAN_TOP')
seisandb = 'DSNC_'
operator = 'GT' # initials of person running program
#----------------------------------------#

wavdbdir = os.path.join(seisantopdir, 'WAV', seisandb)
if os.path.isdir(wavdbdir):
    print(wavdbdir, ' exists.')
else:
    print(wavdbdir, ' not found. Exiting')
    exit()
readbdir = os.path.join(seisantopdir, 'REA', seisandb)
if os.path.isdir(readbdir):
    print(readbdir, ' exists.')
else:
    print(readbdir, ' not found. Did you run makerea? Exiting')
    exit()

# Commands normally entered by keyboard when running autoreg
cwd = os.getcwd()
autoreg_infile = os.path.join(cwd, 'autoreg_commands.txt')
fptr = open(autoreg_infile, 'w')
fptr.write('L\n') # L, R, or D class
fptr.write('m\n') # c (copy) or m (move)
fptr.write('%s\n' % seisandb) # SEISAN_DB you are registering into
fptr.write('%s\n' % operator) # SEISAN_DB you are registering into
fptr.close()

test_mode = True
n_good = 0 # number of good wavfiles processed so far
fbad = open('badwavfiles.txt', 'w') # log bad WAV files in this file

# --- TO LOOP OVER A WAV/DB/YYYY/MM archive ---#
# --- For a single directory, can just call autoreg_this_dir directly ---#
# Loop over year directories
yeardirs = glob.glob(os.path.join(wavdbdir, '[12][0-9][0-9][0-9]'))
for yeardir in sorted(yeardirs):
    print('Processing ',yeardir)
    yyyy = os.path.basename(yeardir)
    # Loop over month directories
    monthdirs = glob.glob(os.path.join(yeardir, '[01][0-9]'))
    for monthdir in sorted(monthdirs):
        print('Processing ',monthdir)
        autoreg_this_dir(monthdir, yyyy, autoreg_infile, fbad, n_good)

fbad.close()
