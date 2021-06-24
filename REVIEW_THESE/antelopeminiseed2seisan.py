from obspy import read
#from obspy.core import Stream
import os.path
import numpy as np
import datetime
import sys
#import matplotlib.pyplot as plt
import struct
import os
import glob
#import time
#import re
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
#sys.path.insert(0,'~/scripts')
#import libMVOarchive

def process(cwd, seedfiles, seisanfile):
    for seedfile in seedfiles:
        print("Reading %s" % seedfile)
        stnew = read(seedfile)
        #stnew = read(cwd + "/" + seedfile)
        print(stnew)
        try:
            st
        except:
            st = stnew
        else:
            st += stnew
        #stnew.plot(outfile="%s_%s_%s.png" % (seisanfile,stnew[0].stats.station,stnew[0].stats.channel))
    st.plot(outfile="%s.png" % seisanfile)
    st.write(seisanfile, format="MSEED")
    return

def main(argv):
    # Usage: PROGRAM sourcedir seisan_top seisan_dbname startyear startday endyear endday

    # Defaults
    sourcedir = "/Users/tom/Desktop/IRIG" # Where the miniseed files sit under year/julday directories
    sourcedir = "/Volumes/data/Montserrat/MontserratSeismicMastering/seed/IRIG_timecorrected"
    #sourcedir = "/raid/data/Montserrat/MontserratSeismicMastering/seed/IRIG_timecorrected"
    seisan_top = "/Users/tom/Desktop/seismo" # Where Seisan is installed
    seisan_top = "/Users/thompsong/seismo"
    #seisan_top = "/home/t/thompsong/seismo"
    seisan_dbname = "ASNE_" # The Seisan database name
    startyear = "1996"
    startday = "1"
    endyear = "2004"
    endday = "366"
    endyear = "1996"
    endday = "1"

    # Process command line arguments
    if len(argv) > 1:
        sourcedir = argv[1]
        if sourcedir == '.':
            sourcedir = os.getcwd()
    if len(argv) > 2:
        seisan_top = argv[2]
        if seisan_top == '.':
            seisan_top = os.getcwd()
    if len(argv) > 3:
        seisan_dbname = argv[3]
    if len(argv) > 4:
        startyear = argv[4]
    if len(argv) > 5:
        startday = argv[5]
    if len(argv) > 6:
        endyear = argv[6]
        
    if len(argv) > 7:
        endday = argv[7]
        
    startyearjdaystr = startyear + startday.zfill(3)
    endyearjdaystr = endyear + endday.zfill(3)

    # check source directory exists
    yearsourcedir = sourcedir + "/" + startyear
    if os.path.isdir(yearsourcedir):
        print("Will process files from %s" % yearsourcedir)
    else:
        print("Source directory (%s) not found. Exiting." % yearsourcedir)
        return

    # make sure database directory exists in REA and WAV    
    seisan_rea = seisan_top + "/" + "REA" + "/" + seisan_dbname
    if not os.path.exists(seisan_rea):
        reapath = Path(seisan_rea)
        reapath.mkdir(parents=True)
    seisan_wav = seisan_top + "/" + "WAV" + "/" + seisan_dbname
    if not os.path.exists(seisan_wav):
        wavpath = Path(seisan_wav)
        wavpath.mkdir(parents=True)            
    
    # pattern to match
    pattern = '*.mseed'

    # Find & process the files
    if os.path.isdir(sourcedir):
        # directory
        os.chdir(sourcedir)
        # see if it has year subdirectories
        years = glob.glob('[0-9]' * 4)
        if len(years) > 0:
            os.chdir(sourcedir)
            # loop over years
            years.sort()
            for year in years:
                if(year<startyear):
                    continue
                if(year>endyear):
                    return
                os.chdir(sourcedir + "/" + year)
                days = glob.glob('[0-9]' * 3)
                days.sort()
                for day in days:
                    #print(year)
                    #print(day)
                    yearjdaystr = year + day.zfill(3)
                    if(yearjdaystr >= startyearjdaystr and yearjdaystr <= endyearjdaystr ):
                        thisdate = datetime.datetime.strptime(yearjdaystr, '%Y%j').date()
                        jdaystr = thisdate.strftime('%d')
                        print("Processing %s" % jdaystr)
                        month = thisdate.strftime('%m')
                        targetdir = seisan_wav + "/" + year + "/" + month
                        print("Saving to %s" % targetdir)
                        os.chdir(sourcedir + "/" + year + "/" + day)
                        cwd = os.getcwd()
                        files = glob.glob(pattern)
                        files.sort()
                        for file in files:
                            #print(file)
                            parts = file.split(".")
                            #print(parts[2])
                            eventpattern = "*%s.mseed" % parts[2]
                            sameeventfiles = glob.glob(eventpattern)
                            sameeventfiles.sort()
                            #print(sameeventfiles)
                            datetimeparts = parts[2].split("_")
                            hh=datetimeparts[2]
                            mi=datetimeparts[3]
                            ss=datetimeparts[4]
                            seisan_wav_yearmonth = seisan_top + "/" + "WAV" + "/" + seisan_dbname + "/" + year + "/" + month
                            seisan_rea_yearmonth = seisan_top + "/" + "REA" + "/" + seisan_dbname + "/" + year + "/" + month
                            print("Seisan WAV dir: %s" % seisan_wav_yearmonth)
                            numchans = len(sameeventfiles)
                            seisanfilename = seisan_wav_yearmonth + "/" + "%s-%s-%s-%s%s-%sS.%s_%03d" % (year,month,thisdate.strftime('%d'),hh,mi,ss,seisan_dbname,numchans)
                            if not os.path.exists(seisanfilename):
                                print(" ")
                                if not os.path.exists(seisan_wav_yearmonth):
                                    swympath = Path(seisan_wav_yearmonth)
                                    swympath.mkdir(parents=True)
                                if not os.path.exists(seisan_rea_yearmonth):
                                    srympath = Path(seisan_rea_yearmonth)
                                    srympath.mkdir(parents=True)
                                print("Seisan filename: %s" % seisanfilename)
                                for sameeventfile in sameeventfiles:
                                    print("    - %s" % sameeventfile)
                                process(cwd, sameeventfiles, seisanfilename)
                        os.chdir('..')
            os.chdir('..')
        else: # no year directories, just glob local directory
            files = glob.glob(pattern)
            for file in files:
                process(file, targetdir)
    else:
        file = sourcedir
        if os.path.exists(file):
            process(file, targetdir)

    
if __name__ == "__main__":
    # Usage PROGRAMNAME sourcedir seisan_top seisan_dbname
    main(sys.argv)
