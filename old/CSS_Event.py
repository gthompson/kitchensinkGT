#! /usr/bin/env python
""" CSS Event class, Glenn Thompson 2010/05/01
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
class CSS_Event(object):
    NULLSTRING = '-'
    NULLID = -1
    NULLTIME = -9999999999.99900
    NULLREAL = -999.0000
    NULLLOC = '-'
    NULLMAG = -99.99
    NULLM = -999.0

    def __init__(self, evid=1, evname=None, prefor=None, auth=None, commid=None, lddate=None):
        """ initialise an empty Event object """
	print "Creating CSS Event object"
        self.evid = evid
	self.evname = evname
        self.prefor = prefor 
	self.auth = auth
        self.commid = commid
	self.lddate = lddate
	self.origins = list()

    def __str__(self):
        """ print out the attributes of an event object """
	str = "CSS Event object with attributes:\n"
	str += "\tevid: %d\n" % self.evid
	if self.evname:
		str += "\tevname: %s\n" % self.evname
	if self.prefor:
		str += "\tprefor: %d\n" % self.prefor
	if self.auth:
		str += "\tauth: %s\n" % self.auth
	if self.commid:
		str += "\tcommid: %d\n" % self.commid
	if self.lddate:
		str += "\tlddate: %f\n" % self.lddate
	return str

    def add_origin(self, originobject):
        self.origins.append(originobject)

#    def read_event(dbname, evid):

     def parse(self, line):
    	"""
	Parse a line of a CSS Event table and return a CSS_Event object

	The format of an Antelope/CSS3.0 schema event table is taken from $ANTELOPE/data/schemas/css3.0
	(part of the commercial Antelope system)

        Fields ( evid evname prefor auth commid lddate )
		evid [0-7] - integer(%8ld), null -1, evid>0
		evname [8-22] - string(%-15s), null -
		prefor [23-30]  - integer(%8ld), null -1, prefor>0
		auth [31-45] - string(%-15s), null -
		commid [46-53] - integer(%8ld), null -1, commid>0
		lddate [54-70] - time(%17.5f), null -9999999999.99900, epoch seconds
    	"""

	try:
		self.evid = int(line[0:8])
		self.evname = line[8:23].strip().decode()
		self.prefor = int(line[23:31])
		self.auth = line[31:46].strip().decode()
		self.commid = int(line[46:54])
		self.lddate = float(line[54:])

		# check for nulls
		if self.evid == NULLID:
			self.evid = None
		if self.evname == NULLSTRING:
			self.evname = None
		if self.prefor == NULLID
			self.prefor = None
		if self.auth == NULLSTRING:
			self.auth = None
		if self.commid == NULLID
			self.commid = None
		if self.lddate == NULLTIME
			self.lddate = None

	except:
		print "Error reading event from line: %s" % line



if __name__ == "__main__":
    evid = 1
    e = CSS_Event()
    print e
    e = CSS_Event(prefor = 2, auth='Glenn')
    print e

