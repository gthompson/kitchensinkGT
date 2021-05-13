#!/usr/bin/env python
# coding: utf-8

# Glenn Thompson, Feb 2021

# imports
import pandas as pd
import os
import glob

STATSDIR = "eventStats"
allDF = pd.DataFrame()
allEventFiles = sorted(glob.glob(os.path.join(STATSDIR, 'Event.2*.stats.csv') ))
numEvents = len(allEventFiles)
count = 0
for thisEventFile in allEventFiles:
    count = count + 1
    print("Processing %s (%d of %d)" % (thisEventFile, count, numEvents) )
    thisDF = pd.read_csv(thisEventFile)
    if not allDF.empty:
        allDF = allDF.append(thisDF)
    else:
        allDF = thisDF

allDF.to_csv(os.path.join(STATSDIR, 'AllEvents.stats.csv' ))

