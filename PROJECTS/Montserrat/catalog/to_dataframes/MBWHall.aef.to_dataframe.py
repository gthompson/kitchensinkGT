#!/usr/bin/env python
import pandas as pd
df = pd.read_csv('MBWHall.aef.wronglinelengthsremoved.txt', sep='\s+')
print(df)
df.columns = ['date', 'time', 'subclass', 'amp', 'eng', 'peakf', 'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'duration']
print(df)
df.to_csv('MBWH.aef.fixed.csv')
df.to_pickle('MBWH.aef.fixed.pkl')

