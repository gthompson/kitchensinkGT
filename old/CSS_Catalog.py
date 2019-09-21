#! /usr/bin/env python
""" CSS Catalog class, Glenn Thompson 2010/05/01
A Python class & reader for Antelope CSS event tables without using Antelope/Datascope 
Glenn Thompson 2014-02-25

To do:
- Goal for now is to load event and origin tables into a CSSCatalog object
(which I have yet to define)
- Another function is needed to link event and origin tables
- Then another is needed to convert CSSCatalog into ObsPy Catalog
- Then I will want to add Eventrate class, and methods from my MATLAB
catalog/event/eventrate classes
"""

import os
import numpy as np
from obspy import UTCDateTime
class CSS_Catalog(object):
    nullstring = '-'
    nullid = -1
    nulltime = -9999999999.99900
    nullreal = -999.0000
    nullloc = '-'
    nullmag = -99.99
    nullmb = -999.0

    def __init__(self):
        """ initialise an empty Event object """
	print("Creating CSS Catalog object")
	self.events = list()
	self.source = None

    def __str__(self):
        """ print out the attributes of a catalog object """
	str = "CSS Catalog object:\n"
	str += "\tnumber of events: %d\n" % len(self.events)
	if self.source:
		str += "\tsource: %s\n" % self.source
	return str

    def read(self, dbname):
    	"""
    	Read an Antelope CSS3.0 event database table and return a CSS_Event object.

	The format of an Antelope/CSS3.0 schema event table is taken from $ANTELOPE/data/schemas/css3.0
	(part of the commercial Antelope system)
	"""

    	# initialize
    	c = CSS_Catalog()
    
    	# check if event table exists
    	event_table = os.path.join(dbname + '.event')
    	if(os.path.exists(event_table)):
		read_event_table(self, event_table)
	else:
		print("Event table %s does not exist" % event_table)  
    		origin_table = os.path.join(dbname + '.origin')
		read_origin_table(self, origin_table)


    def read_table(table):
	lines = list()
    	with open(table, "r") as fh:
		print("Opening %s" % table)
        	lines = fh.readlines()
		print("- got %d lines" % len(lines))
	return lines

    def read_event_table(self, table):
	lines = read_table(table)
	line_num = 0
    	for line in lines:
		try:
			line_num += 1
			print("processing line %d" % line_num)
			e = CSS_Event()
			e.parse(line)
			self.events.append(e)
		except:
			print("Error reading event at line %d" % line_num)

    def read_origin_table(self, table):
	lines = read_table(table)
	line_num = 0
    	for line in lines:
		try:
			line_num += 1
			print("processing line %d" % line_num)
			o = CSS_Origin()
			o.parse(line)
			self.origins.append(o)
		except:
			print("Error reading origin at line %d" % line_num)


if __name__ == "__main__":
    #c = CSS_Catalog()
    #print c 
    o = CSS_Origin()
    print(o)
    e = CSS_Event()
    print(e)
    c.read('/opt/antelope/data/db/demo/demo')
    print(c)

