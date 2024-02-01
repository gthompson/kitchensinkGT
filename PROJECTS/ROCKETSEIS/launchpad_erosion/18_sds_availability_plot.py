#!/usr/bin/env python
# coding: utf-8

# # Plot SDS data availability

import sys
if len(sys.argv) != 2:
    print(f"{sys.argv[0]} SDS_DIR")
else:
    SDS_DIR = sys.argv[1]
import header
paths = header.setup_environment()
import os
import pandas as pd
from obspy import UTCDateTime
import SDS
import libWellData as LLE

transducersDF = LLE.get_transducers_dataframe(paths)
startover = True

# Availability by day
startdate = UTCDateTime(2022,7,21)
enddate = UTCDateTime(2022,12,3)
#sdsobj = SDS.SDSobj(paths['SDS_TOP'])
sdsobj = SDS.SDSobj(SDS_DIR)
sdsbase = os.path.basename(SDS_DIR)
trace_ids = None
availabilityDIR = os.path.join(paths['outdir'], 'availability')
availabilityCSV = os.path.join(availabilityDIR, f"{sdsbase}.csv")
wellPNG=os.path.join(availabilityDIR,f"{sdsbase}_well.png")
SApng=os.path.join(availabilityDIR,f"{sdsbase}_SA.png")
if startover:
    if os.path.isfile(availabilityCSV):
        os.remove(availabilityCSV)
if os.path.isfile(wellPNG):
    os.remove(wellPNG)
if os.path.isfile(SApng):
    os.remove(SApng)

if os.path.isfile(availabilityCSV): # would have to change startover to False to get here
    availabilityDF = pd.read_csv(availabilityCSV, index_col=None)
    trace_ids = availabilityDF.columns[1:]
else:
    a = sdsobj._sds_percent_available_per_day(startdate, enddate, trace_ids=trace_ids, speed=1)
    print(a)
    availabilityDF, trace_ids = sdsobj._sds_percent_available_per_day(startdate, enddate, trace_ids=trace_ids, speed=1)
    availabilityDF.to_csv(availabilityCSV, index=False)

def reorder_trace_ids(df, ordered_ids):
    print(ordered_ids)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    d = df['date']
    df = df.loc[:, ordered_ids]
    df.insert(0, 'date', d)
    return df

well_availabilityDF = availabilityDF.copy()
for id in trace_ids:
    if id[0] !='6':
        well_availabilityDF.drop(labels=id, axis=1, inplace=True)
well_availabilityDF = reorder_trace_ids(well_availabilityDF, transducersDF['id'].to_list())
sdsobj.plot_availability(well_availabilityDF, outfile=wellPNG, labels=transducersDF.serial.to_list())

SA_availabilityDF = availabilityDF.copy()
for id in trace_ids:
    if id[0] =='6':
        SA_availabilityDF.drop(labels=id, axis=1, inplace=True)
sdsobj.plot_availability(SA_availabilityDF, outfile=SApng)


# Falcon 9 Block 5 | Starlink Group 4-35		September 24, 2022	SLC-40	23:32:10 UTC
