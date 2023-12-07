#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 18:27:02 2021

@author: thompsong
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime as dt
import glob
dropboxpath = '/Users/thompsong/Dropbox/PROFESSIONAL/RESEARCH/3 Project Documents/201602 Rocket Seismology/db/wfdata/rocketevents'

os.chdir(dropboxpath)
eventdirs = sorted(glob.glob("20*"))
df = pd.DataFrame()
for eventdir in eventdirs:
    tryfile = os.path.join(eventdir, 'summary_seismic_3c.csv')
    if os.path.exists(tryfile):
        print('Loading ', tryfile)
        thisdf = pd.read_csv(tryfile)
        print(thisdf)
        df.append(thisdf)
    
print(df)
yyyymmdd = []
for c in range(len(df)):
    yyyy = df.iloc[c]['eventdir'][0:4]
    mm = df.iloc[c]['eventdir'][4:6]
    dd = df.iloc[c]['eventdir'][6:8]
    this_yyyymmdd = dt(int(yyyy), int(mm), int(dd))
    yyyymmdd.append(this_yyyymmdd)
print(yyyymmdd)

#
df.insert(0, 'date', yyyymmdd)
print(df)
#df.plot()

df.columns = df.columns.str.replace(' ', '')
print(df.columns)


df.plot(x="date", y="PGD", kind='scatter')
plt.ylabel('Peak Ground Displacement (m)')

df.plot(x="date", y="PGV", kind='scatter')
plt.ylabel('Peak Ground Velocity (m/s)')

df.plot(x="date", y="PGA", kind='scatter')
plt.ylabel('Peak Ground Acceleration (m/s^2)')