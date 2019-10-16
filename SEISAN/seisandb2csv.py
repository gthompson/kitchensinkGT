#!/Users/thompsong/miniconda3/bin/python
# /Users/thompsong/src/volcanoObsPy/SEISAN/seisandb2csv.py
# Traverses a Seisan database structure by REA/year/month finding all the S-files and
# corresponding WAV files, and generates a summary line for each trace.
# Will generate nothing for an Sfile for which no corresponding WAVfile found
# Output files are like DBNAMEcatalogYYYYMM.csv
# These should be readable into Python dataframes
# They can also easily be manipulated with grep and awk, e.g.
# 1. Only print 4th column:
#    awk -F',' '{print $4}' MVOE_catalog*.csv
# 2. Only show 5th column for station MBWH
#    grep MWBH MVOE_catalog*.csv | awk -F',' '{print $5}'
import os, sys, glob, obspy.core
#sys.path.insert(1, '/Users/thompsong/src/volcanoObsPy')
sys.path.insert(1, '/Users/thompsong/src/volcanoObsPy/SEISAN')
import Seisan_Catalog as SC

def generate_monthly_csv(mm):
    sfileslist = sorted(glob.glob(os.path.join(mm, "*")))
    outfile = mm[-13:-8] + 'catalog' + mm[-7:-3] + mm[-2:] + '.csv'
    print("...Creating %s" % outfile)
    fptr = open(outfile,'w+')
    fptr.write('datetime, mainclass, subclass, duration, wavfilepath, sampling_rate, npts, traceNum, traceID, sfilepath\n') 
    for thissfile in sfileslist:
        print(thissfile)
        s = SC.Sfile(thissfile)
        #print(s)
        numwavfiles = len(s.wavfiles)
        if not s.subclass:
            s.subclass = "_"
        if numwavfiles > 0:
            for thiswavfile in s.wavfiles:
                if not os.path.exists(thiswavfile.path):
                    continue
                st = obspy.read(thiswavfile.path)
                tracenum = 0
                for tr in st:
                     duration = tr.stats.npts / tr.stats.sampling_rate
                     fptr.write("%s,%s,%s,%8.3f,%s,%3d,%6d,%02d,%s,%s\n" % (s.filetime, s.mainclass, \
                         s.subclass, duration, \
                         thiswavfile.path, \
                         tr.stats.sampling_rate, tr.stats.npts, \
                         tracenum, tr.id, thissfile ))
                     tracenum += 1
    fptr.close()


print("Type the path to your Seisan directory [default: ./seismo]")
SEISAN_TOP = input(" ? ")
if not SEISAN_TOP:
    SEISAN_TOP = "./seismo"
print("Type the name of your Seisan database [default: MVOE_]")
SEISAN_DB = input(" ? ")
if not SEISAN_DB:
    SEISAN_DB = "MVOE_"

yyyylist = sorted(glob.glob(os.path.join(SEISAN_TOP, 'REA', SEISAN_DB, '????')))
for yyyy in yyyylist:
    mmlist = sorted(glob.glob(os.path.join(yyyy, '??')))
    for mm in mmlist:
        print("Processing %s:" % mm)
        generate_monthly_csv(mm)

