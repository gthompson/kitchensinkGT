#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Convert OVSM data from 2-minute long Guralp Compressed Format files in a YYYY/MM/ directory tree to a SeisComP Data Structure (SDS)
# SDS directory tree structure / file naming convention is one MiniSEED file per SEED ID per day
# https://www.seiscomp.de/seiscomp3/doc/applications/slarchive/SDS.html
# 
# Source data:
#     from Google Drive folder https://drive.google.com/drive/folders/1DWcdFgYl-nmuNAi4TWr4Vqexv6xfcFEh
#     20240523-Datos_Ruiz_2012/GCF/2012 from Carlos Cardona via an email from Yuly. Each month 04.zip, 05.zip, 06.zip is a zip file 
#     There are also SUDS_multiplexados_1 and _2 folders, which are also monthly zip files
#     from OVSM (Manizales observatory)
#
# Author: Glenn Thompson 2024/11/15
#
import os
import glob
import obspy

# Input data extracted to this directory tree structure
YYYY = '2012'
gcfdir = f'/data/OVSM/GCF/{YYYY}'

# network suggested by Felix
network = 'NR'

# Top of SDS 
SDSdir = '/data/OVSM/SDS'

for monthdir in sorted(glob.glob(os.path.join(gcfdir, '??'))):
    for stationdir in sorted(glob.glob(os.path.join(monthdir, '*'))):
        if os.path.isdir(stationdir):
            dayst = obspy.Stream()
            lastday = None
            for filepath in sorted(glob.glob(os.path.join(stationdir, '*.gcf'))):
                if os.path.isfile(filepath):
                    try:
                        st = obspy.read(filepath, format='GCF')
                    except Exception as e:
                        print('\nUnknown format for ',filepath)
                    else:
                        startt = st[0].stats.starttime
                        currentday = startt.day 
                        if currentday != lastday:
                            if len(dayst)>0:
                                for tr in dayst:
                                    sdsfulldir = os.path.join(SDSdir, network, startt.strftime('%Y'), tr.stats.station, tr.stats.channel+'.D')
                                    sdsfullpath = os.path.join(sdsfulldir, f"{tr.id}.D.{startt.strftime('%Y.%j')}")
                                    if not os.path.isdir(sdsfulldir):
                                        os.makedirs(sdsfulldir)
                                    tr.write(sdsfullpath, format='MSEED')
                            dayst = st.copy()
                            for tr in dayst:
                                tr.stats.network = network
                            lastday = currentday
                        else:
                            for tr in st:
                                tr.stats.network = network
                                dayst.append(tr.copy())
                            dayst.merge(method=0, fill_value=0)

            
            


# In[ ]:




