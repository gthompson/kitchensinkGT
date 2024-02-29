#!/usr/bin/env python
import pandas as pd
df = pd.read_csv('mbwh_subclass_magnitude_sorted.dat',sep='\s+')
df.columns = ['year', 'month', 'day', 'hour', 'minute', 'second', 'subclass', 'emag']
df.to_csv('mbwh_subclass_magnitude_sorted.csv')
df.to_pickle('mbwh_subclass_magnitude_sorted.pkl')
