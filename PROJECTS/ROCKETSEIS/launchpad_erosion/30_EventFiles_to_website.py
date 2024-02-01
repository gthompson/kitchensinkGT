#!/usr/bin/env python
# coding: utf-8

# # Segment SDS Archive based on launch event times and generate event web browser
# 
# Launch times come from 'PilotStudy_KSC_Rocket_Launches.xlsx'
# SDS archives are at SDS_TOP and contain data from wells 6I and 6S and from seismo-acoustic stations
# Segmented event waveform files are saved as MiniSEED to EVENT_WAVEFORMS
# 
# 
import header
paths = header.setup_environment()
import os
#import sys
import glob
#import numpy as np
import pandas as pd
from obspy.core import UTCDateTime, read
#import FDSNtools
#import wrappers
#import SDS
import libWellData as LLE

HTML_DIR = '/var/www/html/thompsong/KSC_EROSION/EVENTS'
PNG_DIR = os.path.join(HTML_DIR, 'images')
EVENT_WAVEFORMS = os.path.join(paths['outdir'], 'EVENTS') # must exist, and Excel file must be here
csv_launches_detected = os.path.join(paths['outdir'], 'PilotStudy_KSC_Rocket_Launches_detected.csv')

launchesDF = LLE.removed_unnamed_columns(pd.read_csv(csv_launches_detected, index_col=None))

def group_streams_for_plotting(st):
    groups = {}
    stationsWELL = ['6S', '6I']
    for station in stationsWELL:
        stationStream = st.select(network=station)
        #stationIDS = list(set([tr.id for tr in stationStream]))
        groups[station] = stationStream
    streamSA = st.select(network='FL')
    stationsSA = list(set([tr.stats.station for tr in streamSA]))
    for station in stationsSA:
        stationStream = streamSA.select(station=station)
        #stationIDS = list(set([tr.id for tr in stationStream]))
        groups[station] = stationStream
    #print(groups)
    return groups  

for thisdir in [HTML_DIR, PNG_DIR]:
    if not os.path.isdir(thisdir):
        os.makedirs(thisdir)

def make_event_html(i, row, groups, ext='short',peakmeas=None,units=None):
    stations = groups.keys()
    launchTime = UTCDateTime(row['Date'])    
    lts = launchTime.strftime('%Y%m%dT%H%M%S')
    lts_human = launchTime.strftime('%Y-%m-%d %H:%M:%S')
    htmlfile = os.path.join(HTML_DIR, 'launch_%s_%s.html' % (lts, ext))
    nl = '\n'
    print(f"Writing {htmlfile}")
    contents = """
<html>
<head>
<title>"""
    contents = f"Event {lts_human}</title>\n</head>\n\n<body>{nl}"

    # EVENT INFO
    contents += "<table border=1>\n"
    contents += f"<tr> <th>Event Number:</th> <td>{i}</td> </tr>\n"
    contents += f"<tr> <th>Date/Time:</th> <td>{lts_human} {launchTime.timestamp,} ({launchTime.strftime('%j')}) </td> </tr>{nl}"
    contents += f"<tr> <th>Detection time:</th> <td>{row['detection_time']}</td> </tr>{nl}"
    contents += f"<tr> <th>Rocket/Payload:</th> <td>{row['Rocket_Payload']}</td> </tr>{nl}"
    contents += f"<tr> <th>Launchpad:</th> <td>{row['SLC']}</td> </tr>{nl}"
    contents += f"<tr> <th>Notes:</th> <td>{row['Notes']}</td> </tr>{nl}"
    contents += "</table>\n"
    
    # PLOTS
    contents += "<table border=1>\n"
    for stationname in groups.keys():
        ids = [tr.id for tr in groups[stationname]] 
        pngfile = os.path.join(os.path.basename(PNG_DIR), '%s_%s_%s.png' % (lts, stationname, ext))
        contents += f"<tr> <td><h1>{stationname}</h1></td> <td><a href='{pngfile}'><img src='{pngfile}'></a><td/> {nl}"
        if peakmeas:
            contents += "<td>"
            for id in ids:
                try:
                    thispeakstr = '%.2e' % peakmeas[id]
                    contents += f"{id}: {thispeakstr} {units[id]} <br/>{nl}"    
                except:
                    pass  
            contents += "</td>"
            contents += "</tr> \n"

    contents += "</table>\n"
    
    contents += "</body>\n</html>"
    fptr = open(htmlfile, "w")
    fptr.write(contents)
    fptr.close()
    #print(contents)
    return os.path.basename(htmlfile)

def make_index_html(launchesDF):
    contents0 = """
<html>
<head>
<title>Events</title>
</head>
<body>
<table border=1>
<tr><th>Date</th><th>Rocket/Payload</th><th>Launchpad</th><th>Raw</th><th>Short</th><th>Long</th><th>Strength (um/s)</th></tr>
    """    
    for i, row in launchesDF.iterrows():

        print(i, 'of ', len(launchesDF) )

        launchTime = UTCDateTime(row['Date'])    
        lts = launchTime.strftime('%Y%m%dT%H%M%S')
        lts_human = launchTime.strftime('%Y-%m-%d %H:%M:%S')

        contents0 += f"<tr> <td>{lts_human} ({launchTime.timestamp,} {launchTime.strftime('%j')}) </td> "
        contents0 += f"<td>{row['Rocket_Payload']}</td> <td>{row['SLC']}</td> "
        rawpng = os.path.join(EVENT_WAVEFORMS, row['rawfile'].replace('.pkl','.png'))
        if os.path.isfile(rawpng):
            contents0 += f"<td><a href={rawpng}>raw</a></td> "
        else:
            contents0 += "<td>None</td> "

        HHZpeak = 0 
        for ext in ['short', 'long']:
            if ext=='short':
                pklfile = row['short_file']
            else:    
                pklfile = row['corrected_file']
            if not pklfile:
                continue
            pklfile =  os.path.join(EVENT_WAVEFORMS, pklfile)
            print('pklfile=',pklfile) 
            if os.path.isfile(pklfile):
                st = read(pklfile, 'PICKLE')
                groups = group_streams_for_plotting(st)
                peakmeas = dict()
                units = dict()
                if ext=='short':
                    for tr in st:
                        peakmeas[tr.id] = max(abs(tr.data)) 
                        units[tr.id] = tr.stats.units
                        if tr.stats.channel == 'HHZ':
                            if peakmeas[tr.id] > HHZpeak:
                                HHZpeak = peakmeas[tr.id]

                htmlfile = make_event_html(i, row, groups, ext, peakmeas, units)
                contents0 += f"<td><a href={htmlfile}>{ext}</a></td> "
            else:
                contents0 += f"<td>None</td> "
        HHZpeakstr = '%.2e' % HHZpeak
        contents0 += f"<td>{HHZpeakstr}</td> </tr>"
        contents0 += '\n'
    contents0 += "\n</table>\n</body>\n</html>"
    indexfile = os.path.join(HTML_DIR, 'index.html')
    fptr0 = open(indexfile, "w")
    fptr0.write(contents0)
    fptr0.close()    


#make_html(UTCDateTime('2022-11-03 05:22:00'), ['S39A1', 'BCHH2'])
#print(launchesDF)
make_index_html(launchesDF)

