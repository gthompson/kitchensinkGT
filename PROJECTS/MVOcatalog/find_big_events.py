#!/usr/bin/env python
import pandas as pd
import csv
import obspy
import matplotlib.pyplot as plt
import numpy as np
import trace_quality_control as qc 
import os
import glob


# Load the catalogue
df = pd.read_csv('./csvfiles/MVOE_sfiles200701.csv')
nrows, ncolumns = df.shape
#print(nrows)
#print(df.head)
#print(df.wavfile1)
df2 = pd.read_csv('./csvfiles/MVOE_wavfiles200701.csv')
print(df2.head)
for index,row in df.iterrows():
	#print(index, row)
	print(row.wavfile1)
	# duration and wavvfilepath headers are swapped
	index2 = df2["duration"].str.find(row.wavfile1)
	print(index2)
#        wavfilepath = mvocat.iloc[i].wavfilepath
        #print(wavfilepath)
#        if os.path.exists(wavfilepath):
