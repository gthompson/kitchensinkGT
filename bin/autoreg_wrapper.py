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

def autoreg_this_dir(monthdir, yyyy, autoreg_infile, fbad, n_good):
    # Loop over each WAV file
    wavfiles = glob.glob(os.path.join(monthdir, '%s*S.*' % yyyy ))
    print('Processing %s: %d WAV files found' % (monthdir, len(wavfiles)))
    for wavfile in wavfiles:
        os.system('dirf %s' % wavfile) # Run dirf - produces filenr.lis with only 1 event
        try: # run autoreg
            os.system('autoreg < %s' % autoreg_infile)
            print('Registered %s' % wavfile)
            n_good += 1
        except: # it crashed. add to list of bad WAV files
            print('autoreg crashed on %s' % wavfile)
            fbad.write('%s\n' % wavfile)
    print('registered %d WAV files so far' % n_good)


# main program starts here
seisantopdir = os.getenv('SEISAN_TOP')
seisandb = 'DSNC_'
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
autoreg_infile = 'autoreg_commands.txt'
operator = 'GT' # initials of person running program
fptr = open(autoreg_infile, 'w')
fptr.write('L\n') # L, R, or D class
fptr.write('m\n') # c (copy) or m (move)
fptr.write('%s\n' % seisandb) # SEISAN_DB you are registering into
fptr.write('%s\n' % operator) # SEISAN_DB you are registering into
fptr.close()

n_good = 0 # number of good wavfiles processed so far
fbad = open('badwavfiles.txt', 'w') # log bad WAV files in this file

# --- TO LOOP OVER A WAV/DB/YYYY/MM archive ---#
# --- For a single directory, can just call autoreg_this_dir directly ---#
# Loop over year directories
yeardirs = glob.glob(os.path.join(wavdbdir, '[12][0-9][0-9][0-9]'))
for yeardir in yeardirs:
    print('Processing ',yeardir)
    yyyy = os.path.basename(yeardir)
    # Loop over month directories
    monthdirs = glob.glob(os.path.join(yeardir, '[01][0-9]'))
    for monthdir in monthdirs:
        print('Processing ',monthdir)
        autoreg_this_dir(monthdir, yyyy, autoreg_infile, fbad, n_good)

fbad.close()
