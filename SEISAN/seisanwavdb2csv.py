#!/usr/bin/env python
# Traverses a Seisan database structure by WAV/db/year/month finding all the WAV-files and generate a summary line for each trace
#
# Output files are like DBNAMEwavfilesYYYYMM.csv
#
# These will later be read into Python dataframes
import sys, glob, os
import Seisan_Catalog

SEISAN_DATA = os.environ["SEISAN_DATA"]
if not SEISAN_DATA:
    SEISAN_DATA = "./seismo"
#print(sys.argv)
#print(len(sys.argv))
if len(sys.argv)>1:
    SEISAN_DB = sys.argv[1]
else:
    print("Type the path to your Seisan directory [default: %s]" % SEISAN_DATA)
    SEISAN_DATA_2 = input(" ? ")
    if SEISAN_DATA_2:
        SEISAN_DATA = SEISAN_DATA_2
    print("Type the name of your Seisan database [default: MVOE_]")
    SEISAN_DB = input(" ? ")
    if not SEISAN_DB:
        SEISAN_DB = "MVOE_"

yyyylist = sorted(glob.glob(os.path.join(SEISAN_DATA, 'WAV', SEISAN_DB, '????')))
print(yyyylist)
for yyyy in yyyylist:
    mmlist = sorted(glob.glob(os.path.join(yyyy, '??')))
    for mm in mmlist:
        #print("Processing %s:" % mm)
        Seisan_Catalog.generate_monthly_wav_csv(mm)

