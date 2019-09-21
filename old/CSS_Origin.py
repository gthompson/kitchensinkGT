#! /usr/bin/env python
""" CSS Origin class, Glenn Thompson 2010/05/01
A Python class & reader for Antelope CSS origin tables without using Antelope/Datascope 
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
class CSS_Origin(object):
    NULLSTRING = '-'
    NULLID = -1
    NULLTIME = -9999999999.99900
    NULLREAL = -999.0000
    NULLLOC = '-'
    NULLMAG = -99.99
    NULLM = -999.0
    NULLINT = -1

    def __init__(self, lat=None, lon=None, depth=None, time=None, orid=1, evid=None, jdate=None, ndef=None, nass=None, 
		ndp=None, grn=None, srn=None, etype=None, review=None, mbid=None, mb=None, msid=None, ms=None, mlid=None, ml=None, 
		depdp=None, dtype=None, algorithm=None, auth=None, commid=None, lddate=None):
        """ initialise an empty Event object """
	print "Creating CSS Origin object"
	self.lat = lat
	self.lon = lon
	self.depth = depth
	self.time = time
        self.orid = orid
        self.evid = evid
	self.jdate = jdate
	self.ndef = ndef
	self.nass = nass
	self.ndp = ndp
	self.grn = grn
	self.srn = srn
	self.etype = etype
	self.review = review
	self.mbid = mbid
	self.mb = mb
	self.msid = msid
	self.ms = ms
	self.mlid = mlid
	self.ml = ml
        self.depdp = depdp 
        self.dtype = dtype 
        self.algorithm = algorithm 
	self.auth = auth
        self.commid = commid
	self.lddate = lddate
	self.arrivals = list()
	self.netmags = list()
	self.origerr = None

    def __str__(self):
        """ print out the attributes of an origin object """
	str = "CSS Origin object with attributes:\n"

	if self.lat:
		str += "\tlat: %s\n" % self.lat
	if self.lon:
		str += "\tlon: %s\n" % self.lon
	if self.depth:
		str += "\tdepth: %s\n" % self.depth
	if self.time:
		str += "\ttime: %d\n" % self.time
	str += "\torid: %d\n" % self.orid
	if self.evid:
		str += "\tevid: %s\n" % self.evid
	if self.jdate:
		str += "\tjdate: %s\n" % self.jdate
	if self.nass:
		str += "\tnass: %s\n" % self.nass
	if self.ndef:
		str += "\tndef: %s\n" % self.ndef
	if self.ndp:
		str += "\tndp: %s\n" % self.ndp
	if self.grn:
		str += "\tgrn: %s\n" % self.grn
	if self.srn:
		str += "\tsrn: %s\n" % self.srn
	if self.etype:
		str += "\tetype: %s\n" % self.etype
	if self.review:
		str += "\treview: %s\n" % self.review
	if self.mb:
		str += "\tmb: %s\n" % self.mb
	if self.ms:
		str += "\tms: %s\n" % self.ms
	if self.ml:
		str += "\tml: %s\n" % self.ml
	if self.depdp:
		str += "\tdepdp: %s\n" % self.depdp
	if self.dtype:
		str += "\tdtype: %s\n" % self.dtype
	if self.algorithm:
		str += "\talgorithm: %s\n" % self.algorithm
	if self.auth:
		str += "\tauth: %s\n" % self.auth
	if self.commid:
		str += "\tcommid: %d\n" % self.commid
	if self.lddate:
		str += "\tlddate: %f\n" % self.lddate
	if self.arrivals:
		str += "\tarrivals: %f\n" % len(self.arrivals)
	if self.netmags:
		str += "\tnetmags: %f\n" % len(self.netmags)
	return str

    def add_arrival(self, arrivalobject):
        self.arrivals.append(arrivalobject)

    def add_netmag(self, netmagobject):
        self.netmags.append(netmagobject)

    def add_origerr(self, origerrobject):
        self.origerr = origerrobject 

#   def read_event(dbname, evid):

    def parse(self, line):
    	"""
	Parse a line of a CSS Origin table and return a CSS_Origin object

	The format of an Antelope/CSS3.0 schema origin table is taken from $ANTELOPE/data/schemas/css3.0
	(part of the commercial Antelope system)

        Fields ( lat lon depth time orid evid jdate nass ndef ndp grn srn etype review depdp dtype mb mbid ms msid ml mlid algorithm auth commid lddate )
		lat - real(%9.4f), null -999.0000, lat >= -90.0 && lat <= 90.0
		lon - real(%9.4f), null -999.0000, lon >= -180.0 && lon <= 180.0
		depth - real(%9.4f), null -999.0000, depth >= 0.0 && depth < 1000.0, kilometres
		time - time(%17.5f), null -9999999999.99900, epoch seconds
                orid - integer(%8ld), null -1, orid>0
                evid - integer(%8ld), null -1, evid>0
		jdate - yearday(%8ld), null -1, year > 0000, day > 000 (obsolete because of time)
		nass - integer(%4ld), null -1, nass >=0 (number of assoc phases)
		ndef - integer(%4ld), null -1, ndef >=0 (number of locating phases)
		ndp - integer(%4ld), null -1, ndp >=0 (number of depth phases)
		grn - integer(%8ld), null -1, grn > 0 (geographic region number)
		srn - integer(%8ld), null -1, srn > 0 (seismic region number)
		etype - string(%-2s), null -, etype =~ /qb|eq|me|ex|o|l|r|t|f/
		review - string(%-4s), null - (y=inspected, auto, pre, fin)
		depdp - real(%9.4f) null -999.0000,  "depdp >= 0.0 && depdp < 1000.0"
		dtype - string(%-1s), null -, "dtype =~ /f|d|r|g/"
		mb - real(%7.2f), null -999.0000, body wave magnitude
                mbid - integer(%8ld), null -1, mbid>0
		ms - real(%7.2f), null -999.0000, surface wave magnitude
                msid - integer(%8ld), null -1, msid>0
		ml - real(%7.2f), null -999.0000, local magnitude
                mlid - integer(%8ld), null -1, mlid>0
		algorithm - string(%-15s), null - (Processing algorithm used)
		auth - string(%-15s), null -
		commid - integer(%8ld), null -1, commid>0
		lddate - time(%17.5f), null -9999999999.99900, epoch seconds
    	"""

	try:
		self.lat = float(line[0:9])
		self.lon = float(line[9:18])
		self.depth = float(line[18:27])
		self.time = float(line[27:44])
		self.orid = int(line[44:52])
		self.evid = int(line[52:60])
		self.jdate = int(line[60:68])
		self.nass = int(line[68:72])
		self.ndef = int(line[72:76])
		self.ndp = int(line[76:80])
		self.grn = int(line[80:88])
		self.srn = int(line[88:96])
		self.etype = line[96:98].strip().decode()
		self.review = line[98:102].strip().decode()
		self.depdp = float(line[102:111])
		self.dtype = line[102:103].strip().decode()
		self.mb = float(line[103:110])
		self.mbid = int(line[110:118])
		self.ms = float(line[118:125])
		self.msid = int(line[125:133])
		self.ml = float(line[133:140])
		self.mlid = int(line[140:148])
		self.algorithm = line[148:163].strip().decode()
		self.auth = line[163:178].strip().decode()
		self.commid = int(line[178:186])
		self.lddate = float(line[186:203])

		# check for nulls
		if self.lat == NULLREAL:
			self.lat = None
		if self.lon == NULLREAL:
			self.lon = None
		if self.depth == NULLREAL:
			self.depth = None
		if self.time == NULLTIME:
			self.time = None
		if self.orid == NULLID:
			self.orid = None
		if self.nass == NULLINT:
			self.nass = None
		if self.ndef == NULLINT:
			self.ndef = None
		if self.ndp == NULLINT:
			self.ndp = None
		if self.grn == NULLINT:
			self.grn = None
		if self.srn == NULLINT:
			self.srn = None
		if self.etype == NULLSTRING:
			self.etype = None
		if self.review == NULLSTRING:
			self.review = None
		if self.depdp == NULLREAL:
			self.depdp = None
		if self.dtype == NULLREAL:
			self.dtype = None
		if self.mb == NULLM:
			self.mb = None
		if self.mbid == NULLID:
			self.mbid = None
		if self.ms == NULLM:
			self.ms = None
		if self.msid == NULLID:
			self.msid = None
		if self.ml == NULLM:
			self.ml = None
		if self.mlid == NULLID:
			self.mlid = None
		if self.algorithm == NULLSTRING:
			self.algorithm = None
		if self.auth == NULLSTRING:
			self.auth = None
		if self.commid == NULLID:
			self.commid = None
		if self.lddate == NULLTIME:
			self.lddate = None

	except:
		print "Error reading event from line: %s" % line



if __name__ == "__main__":
    o = CSS_Origin(orid=5, evid=4, lat=49.2, lon=3.4, depth=65.2, algorithm='manual', ml=3.7, mb=4.0, auth='glenn')
    print o


    """
    Relation netmag
        Fields ( magid net orid evid magtype nsta magnitude uncertainty auth commid lddate )
    
    Relation assoc
        Fields ( arid orid sta phase belief delta seaz esaz timeres timedef azres azdef slores slodef emares wgt vmodel commid lddate )

    Relation arrival
        Fields ( sta time arid jdate stassid chanid chan iphase stype deltim azimuth delaz slow delslo ema rect amp per logat clip fm snr qual auth commid lddate )

    Relation stamag
        Fields ( magid sta arid orid evid phase magtype magnitude uncertainty auth commid lddate )

    Relation origerr
        Fields ( orid sxx syy szz stt sxy sxz syz stx sty stz sdobs smajax sminax strike sdepth stime conf commid lddate )

    """

