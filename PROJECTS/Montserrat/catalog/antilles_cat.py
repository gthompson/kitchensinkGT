#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 17:16:04 2021

@author: thompsong
"""


from obspy import read_events
import os
import basemap
HOME = os.getenv('HOME')
cat = read_events('%s/Dropbox/lesser_antilles_catalog.xml' % HOME)

cat.plot()