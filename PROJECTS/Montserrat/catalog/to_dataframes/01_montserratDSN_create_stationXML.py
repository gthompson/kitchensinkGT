#!/usr/bin/env python
# coding: utf-8

# not for running at MESS2024
# This script is intended to be run on hal9000 at USF seismic lab.
# It parses the MVO station0.HYP file from Seisan to retrieve station coordinates
# It then reads single channel StationXML files created through some unremembered code 2-3 years ago,
# inserts the station coordinaates at the channel and station levels in the Inventory object, and then
# combines all the Inventory objects before saving a StationXML file for the entire MVO DSN.

import os
import sys
import glob
import obspy
from obspy.core.inventory.inventory import Inventory, read_inventory
sys.path.append('..')
import setup_paths
paths = setup_paths.paths
LIB_DIR = os.path.join(paths['PROJECT_DIR'],'src','lib')
sys.path.append(LIB_DIR)
import InventoryTools

stationXMLfile = os.path.join(paths['RESPONSE_DIR'], 'MV.xml')
if os.path.isfile(stationXMLfile):
    exit()

# parse Seisan STATION0.HYP file to get station coordinates
coordinates={}
station0hypfile = '/data/SEISAN_DB/STATION0_MVO.HYP'
fptr= open(station0hypfile , 'r')
lines = fptr.readlines()
for line in lines:
    line.strip()
    if line[2:4]=='MB':
        print(line)
        station=line[2:6].strip()
        #print(line[6:8])
        lat=int(line[6:8]) + int(line[8:10])/60
        if line[10]=='.':
            lat += int(line[11:13])/3600
        else:
            lat += int(line[10:13])/36000
        if line[14]=='S':
            lat=-lat
        lon=int(line[15:17]) + int(line[17:19])/60
        if line[19]=='.':
            lon += int(line[20:22])/3600
        else:
            lon += int(line[19:22])/36000
        if line[22]=='W':
            lon=-lon      
        elev=int(line[23:27])
        #print(f"station={station} lat={lat}, lon={lon}, elev={elev}")
        coordinates[station]={'lat':lat, 'lon':lon, 'elev':elev}
fptr.close()
print(coordinates)

# read stationXML for each trace
stationXMLsrcDIR = '/data/SEISAN_DB/CAL'
stationXMLfiles = sorted(glob.glob(os.path.join(stationXMLsrcDIR, 'station.MV.*..[BSEH]*.xml')))
invAll = None
for xmlfile in stationXMLfiles:
    print(xmlfile)
    this_inv = read_inventory(xmlfile)
    #print(this_inv)
    trace_ids = InventoryTools.inventory2traceid(this_inv, force_location_code='')
    print(trace_ids)
    for seed_id in trace_ids:
        #print(this_inv.get_channel_metadata(id))
        net,sta,chan,loc=seed_id.split('.')
        station_coords = coordinates[sta]

        # - add coordinates
        InventoryTools.modify_inventory(this_inv, seed_id, lat=station_coords['lat'], lon=station_coords['lon'], elev=station_coords['elev'])
        #print(this_inv.get_channel_metadata(id))
    if invAll:
        # combine inventories 
        invAll.extend(this_inv)
    else:
        invAll = this_inv
        
invAll.write(stationXMLfile, format='stationxml')

#invAll.get_channel_metadata('MV.MBBE..BHE')
#trace_ids = InventoryTools.inventory2traceid(invAll, force_location_code='')
#print(trace_ids)

