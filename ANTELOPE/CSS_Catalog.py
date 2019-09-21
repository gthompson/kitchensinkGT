#! /usr/bin/env python
"""
Datascope CSS3.0 Catalog handler
Read databases generated with Antelope and export them to ObsPy Catalog objects
Glenn Thompson, 2014/02/27-2014/03/02.

The Center for Seismic Studies relational database schema version 3.0 is implemented
in Datascope, the relational database management system that is incorporated in Antelope.
These Antelope CSS3.0 databases have tables which describe earthquake catalogs, station
metadata and waveform data.
The earthquake catalog database tables are event, origin, origerr, netmag, stamag, assoc and arrival.

This format is used extensively by PASSCAL campaign projects, and dozens of seismic networks and data centers worldwide. 

The format of an Antelope/CSS3.0 schema database is taken from $ANTELOPE/data/schemas/css3.0
(part of the commercial Antelope system)

In this format an event catalog is described by the tables: event, origin, origerr, netmag, stamag, assoc and arrival.

Currently only event and origin are read in.

Not read in yet are the following:

    netmag
        Fields ( magid net orid evid magtype nsta magnitude uncertainty auth commid lddate )
    
    assoc
        Fields ( arid orid sta phase belief delta seaz esaz timeres timedef azres azdef slores slodef emares wgt vmodel commid lddate )

    arrival
        Fields ( sta time arid jdate stassid chanid chan iphase stype deltim azimuth delaz slow delslo ema rect amp per logat clip fm snr qual auth commid lddate )

    stamag
        Fields ( magid sta arid orid evid phase magtype magnitude uncertainty auth commid lddate )

    origerr
        Fields ( orid sxx syy szz stt sxy sxz syz stx sty stz sdobs smajax sminax strike sdepth stime conf commid lddate )

"""

import os
from obspy import UTCDateTime
import datetime
from obspy.core.event import Catalog, Event, Origin, Magnitude, ResourceIdentifier, CreationInfo
from obspy.geodetics import FlinnEngdahl
VERBOSE = False
ANTELOPE_INSTALLED = False
try:
    from antelope.datascope import *
    ANTELOPE_INSTALLED = True
except:
    # module doesn't exist, deal with it.
    print("antelope.datascope not installed")

class CSS_Catalog(object):
    """
    CSS_Catalog class
    Blueprint for CSS3.0 catalog objects
    """

    def __init__(self):
        """ initialise an empty Event object """
        if VERBOSE:
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

    def import_antelope(self, table, expr=""):
        """ load a catalog using Antelope/Datascope tools """
        """ if an event table is given, it will be joined to origin on prefor """
        """ if an origin table is given, it is assumed there is no event table """
        print("IMPORT ANTELOPE")
        parts = table.split('.')
        dbname = parts[0]
        self.source = dbname
        db = dbopen(dbname, 'r')
        dbo = dbopen_table(db, 'origin')
        if parts[1]=='event':
            dbe = dbopen_table(db, 'event')
            dbj = dbjoin(dbe, dbo)
            dbs = dbsubset(dbj, 'orid == prefor')
        elif parts[1]=='origin':
            dbs = dbo
        nrecs = dbquery(dbs, 'dbRECORD_COUNT')
        for i in range(nrecs):
            dbs[3] = i-1
            e = CSS_Event()
            if parts[1]=='event':
                [e.evid, e.evname, e.prefor, e.auth, e.commid, e.lddate] = \
                 dbgetv(dbs, 'evid', 'evname', 'prefor', 'event.auth', 'event.commid', 
                  'event.lddate')
            elif parts[1]=='origin':
                e.evid = o.evid
                e.prefor = o.orid
            o = CSS_Origin()
            [o.lat, o.lon, o.depth, o.time, o.orid, o.evid, o.jdate, o.nass,
             o.ndef, o.ndp, o.grn, o.srn, o.etype, o.review, o.depdp, 
             o.mb, o.mbid, o.ms, o.msid, o.ml, o.mlid, o.algorithm, 
             o.auth, o.commid, o.lddate] = \
             dbgetv(dbs, 'lat', 'lon', 'depth', 'time', 'orid', 'origin.evid', 'jdate',
              'nass', 'ndef', 'ndp', 'grn', 'srn', 'etype', 'review', 'depdp',
              'mb', 'mbid', 'ms', 'msid', 'ml', 'mlid', 'algorithm',
              'origin.auth', 'origin.commid', 'origin.lddate')
            e.origins.append(o)
def parse_attribute(line, length, ptr, attribute_type):
    """
    function [pend, var] = parse_attribute(line, length, ptr, attribute_type)
    Helper function to parse lines read from a CSS3.0 table
    Inputs:
        line           = line of text (read in from file) to parse
        length         = length of the variable to read
        ptr            = pointer position
        attribute_type = data type, e.g. string, float, int, epoch time

    Outputs:
        pend           = pointer position after reading variable
        var            = attribute value (default:  None)
    """
    # Define Null values used in Datascope
    NULLINT = -1
    NULLREAL = -999.0000
    NULLSTRING = '-'
    NULLTIME = -9999999999.99900
    #NULLMAG = -99.99
    NULLM = -999.0
    var = None                           # The value to return
    pstart = ptr                         # Line position to start reading from
    pend = pstart+length                 # Line position to end reading
    str = line[pstart:pend].strip()      # string variable to hold what is read from line
                                         #        strip it to remove whitespace
    if len(str) > 0:
        if attribute_type=='s' and str!=NULLSTRING: # string type
            var = str
        elif attribute_type=='i':                   # integer type
            try:
                if int(str)!=NULLINT:
                    var = int(str)
            except:
                var = None
        elif attribute_type=='f':                   # float type
            try:
                if float(str)!=NULLREAL:
                    var = float(str)
            except:
                var = None
        elif attribute_type=='t':                   # epoch time type
            try:
                dummy = float(str)
                if dummy!=NULLTIME:
                    var = UTCDateTime(dummy)
            except:
                var = None
        elif attribute_type=='m':                   # magnitude type
            try:
                dummy = float(str)
                if dummy!=NULLM:
                    var = dummy
            except:
                var = None
        else:                                       # unknown type
            var = None
    #print pstart, pend, var
    return pend, var

def read_table(table):
    """
    function lines = read_table(table)
    Helper function to read a CSS3.0 table
    Inputs:
        table      = path of the table to read

    Outputs:
        lines      = lines read from the file
    """
    lines = list()
    with open(table, "r") as fh:
        print("Opening %s" % table)
        lines = fh.readlines()
        print("- got %d lines" % len(lines))
    return lines


class CSS_Catalog(object):
    """
    CSS_Catalog class
    Blueprint for CSS3.0 catalog objects
    """

    def __init__(self):
        """ initialise an empty Event object """
        if VERBOSE:
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

    def import_antelope(self, table, expr=""):
        """ load a catalog using Antelope/Datascope tools """
        """ if an event table is given, it will be joined to origin on prefor """
        """ if an origin table is given, it is assumed there is no event table """
        print("IMPORT ANTELOPE")
        parts = table.split('.')
        dbname = parts[0]
        self.source = dbname
        db = dbopen(dbname, 'r')
        dbo = dbopen_table(db, 'origin')
        if parts[1]=='event':
            dbe = dbopen_table(db, 'event')
            dbj = dbjoin(dbe, dbo)
            dbs = dbsubset(dbj, 'orid == prefor')
        elif parts[1]=='origin':
            dbs = dbo
        nrecs = dbquery(dbs, 'dbRECORD_COUNT')
        for i in range(nrecs):
            dbs[3] = i-1
            e = CSS_Event()
            if parts[1]=='event':
                [e.evid, e.evname, e.prefor, e.auth, e.commid, e.lddate] = \
                 dbgetv(dbs, 'evid', 'evname', 'prefor', 'event.auth', 'event.commid', 
                  'event.lddate')
            elif parts[1]=='origin':
                e.evid = o.evid
                e.prefor = o.orid
            o = CSS_Origin()
            [o.lat, o.lon, o.depth, o.time, o.orid, o.evid, o.jdate, o.nass,
             o.ndef, o.ndp, o.grn, o.srn, o.etype, o.review, o.depdp, 
             o.mb, o.mbid, o.ms, o.msid, o.ml, o.mlid, o.algorithm, 
             o.auth, o.commid, o.lddate] = \
             dbgetv(dbs, 'lat', 'lon', 'depth', 'time', 'orid', 'origin.evid', 'jdate',
              'nass', 'ndef', 'ndp', 'grn', 'srn', 'etype', 'review', 'depdp',
              'mb', 'mbid', 'ms', 'msid', 'ml', 'mlid', 'algorithm',
              'origin.auth', 'origin.commid', 'origin.lddate')
            e.origins.append(o)
            self.events.append(e)
        dbclose(db)

    def read_event_table(self, table):
        """ 
        read a CSS3.0 event table into a catalog object
        """
        # if Antelope is installed, we can use the easier Datascope methods
        #if ANTELOPE_INSTALLED:
        #    import_antelope(table)
        #    return
        lines = read_table(table)
        line_num = 0
        for line in lines:
            line_num += 1
            if VERBOSE:
                print("Processing line %d" % line_num)
            e = CSS_Event()        # instantiate blank event object
            e.parse(line)          # parse line for this event
            self.events.append(e)  # append event to catalog

    def read_origin_table(self, table):
        """ 
        read a CSS3.0 origin table into a catalog object
        only called on a catalog object if no event table is available
        dummy event objects are created to attach the origin objects to
        since origin objects cannot connect directly to catalog objects
        """
        # if Antelope is installed, we can use the easier Datascope methods
        #if ANTELOPE_INSTALLED:
        #    import_antelope(table)
        #    return
        lines = read_table(table)
        line_num = 0
        for line in lines:
            line_num += 1
            if VERBOSE:
                print("processing line %d" % line_num)
            o = CSS_Origin()
            o.parse(line)
            if o.evid is None:
                o.evid = line_num
            e = CSS_Event(evid=o.evid, prefor=o.orid)
            e.origins.append(o)
            self.events.append(e)

    def read(self, dbname):
        """
        Read a CSS3.0 Catalog into a CSS_Catalog object
        """
        # check if event table exists
        # if it does, read it
        # otherwise, read in the origin table instead
        event_table = os.path.join(dbname + '.event')
        if(os.path.exists(event_table)):
            self.read_event_table(event_table)
        else:
            print("Event table %s does not exist" % event_table)  
            origin_table = os.path.join(dbname + '.origin')
            if(os.path.exists(origin_table)):
                self.read_origin_table(origin_table)

    def to_obspy(self):
        """
        Convert CSS_Catalog object into an ObsPy catalog object
        """
        cat = Catalog()
        cat.description = self.source
        for ce in self.events:
            e = Event()
            #e.resource_id = ResourceIdentifier(id="smi:local/event/%d" % ce.evid)
            e.resource_id = ResourceIdentifier(prefix = "event", id="%d" % ce.evid) 
            e.resource_id.convertIDToQuakeMLURI()
            for co in ce.origins:   
                if co.orid == ce.prefor:
                    if co.etype == 'a' or co.etype == 't':
                        e.event_type = 'earthquake'
                        e.event_descriptions.append(EventDescription('volcano-tectonic'))
                    elif co.etype == 'b' or co.etype == 'l':
                        e.event_type = 'earthquake'
                        e.event_descriptions.append(EventDescription('long-period'))
                    elif co.etype == 'h':
                        e.event_type = 'earthquake'
                        e.event_descriptions.append(EventDescription('hybrid'))
                    elif co.etype == 'r':
                        e.event_type = 'debris avalanche'
                        e.event_descriptions.append(EventDescription('rockfall'))
                    elif co.etype == 'c':
                        e.event_type = 'volcanic eruption'
                        e.event_descriptions.append(EventDescription('explosion quake'))
                    elif co.etype == 'e':
                        e.event_type = 'debris avalanche'
                        e.event_descriptions.append(EventDescription('lp-rf'))
                    elif co.etype == 'R':
                        e.event_type = 'earthquake'
                        e.event_descriptions.append(EventDescription('Regional'))
                    elif co.etype == 'L':
                        e.event_type = 'earthquake'
                        e.event_descriptions.append(EventDescription('Local'))
                    elif co.etype == 'T' or co.etype == 'D':
                        e.event_type = 'earthquake'
                        e.event_descriptions.append(EventDescription('Teleseismic'))
                    e.creation_info = CreationInfo(author='co.auth', creation_time=UTCDateTime(datetime.datetime.fromtimestamp(co.lddate)))
                o = Origin()
                #o.resource_id = ResourceIdentifier("%d" % co.orid)
                #o.id = ResourceIdentifier("origin%d" % co.orid)
                #o.id = ResourceIdentifier(id="smi:local/origin/%d" % co.orid)
                o.resource_id = ResourceIdentifier(prefix = "origin", id="%d" % co.orid) 
                o.resource_id.convertIDToQuakeMLURI()
                o.time = co.time
                o.latitude = co.lat
                o.longitude = co.lon
                o.depth = co.depth * 1000
                o.depth_type = "from location"
                if co.review == 'y':
                    o.evaluation_mode = "manual"
                    o.evaluation_status = "final"
                else:
                    o.evaluation_status = "preliminary"
                o.arrivals = list()
                o.origin_type = 'hypocenter'
                o.creation_info = CreationInfo(author='co.auth', creation_time=UTCDateTime(datetime.datetime.fromtimestamp(co.lddate)))
                o.region = FlinnEngdahl().get_region(o.longitude, o.latitude)
                if co.mb:
                    m = Magnitude()
                    m.mag = co.mb
                    m.magnitude_type = "mb"
                    m.origin_id = o.resource_id
                    e.magnitudes.append(m)
                if co.ms:
                    m = Magnitude()
                    m.mag = co.ms
                    m.magnitude_type = "ms"
                    m.origin_id = o.resource_id
                    e.magnitudes.append(m)
                if co.ml:
                    m = Magnitude()
                    m.mag = co.ml
                    m.magnitude_type = "ml"
                    m.origin_id = o.resource_id
                    e.magnitudes.append(m)
                e.origins.append(o)
                cat.append(e)
            
                # also included could be: arrivals, picks, amplitudes, station_magnitudes, 
                #      focal_mechanisms
        return cat

    def to_antelope(self, dbname):
        """ export CSS_Catalog as an Antelope origin table """
        if not ANTELOPE_INSTALLED:
            return
        db = dbopen(dbname, 'a+')
        dbo = dbopen_table(db + '.origin', 'a+')
        for e in self.events:
            # create a new row?
            # use dbaddv or dbputv ?
            # origin table
            for o in e.origins:
                dbputv(dbo, 'orid', o.evid)
                dbputv(dbo, 'time', o.time) # convert datetime to epoch time
                dbputv(dbo, 'lon', o.lon)
                dbputv(dbo, 'lat', o.lat)
                dbputv(dbo, 'depth', o.depth)
                dbputv(dbo, 'ml', o.mag)
            # event table
            dbputv(dbe, 'evid', e.evid)
            dbputv(dbe, 'prefor', e.prefor)
            dbclose(db)

class CSS_Event(object):
    """
    Blueprint for a CSS_Event object

    The format of an Antelope/CSS3.0 schema event table is taken from $ANTELOPE/data/schemas/css3.0
    (part of the commercial Antelope system)

    Fields:
    evid [0-7] - integer(%8ld), null -1, evid>0
    evname [8-22] - string(%-15s), null -
    prefor [23-30]  - integer(%8ld), null -1, prefor>0
    auth [31-45] - string(%-15s), null -
    commid [46-53] - integer(%8ld), null -1, commid>0
    lddate [54-70] - time(%17.5f), null -9999999999.99900, epoch seconds
    """

    def __init__(self, evid=1, evname=None, prefor=None, auth=None, commid=None, lddate=None):
        """ initialise an empty Event object """
        if VERBOSE:
            print("Creating CSS Event object")
        self.evid = evid
        self.evname = evname
        self.prefor = prefor 
        self.auth = auth
        self.commid = commid
        self.lddate = lddate
        self.origins = list()

    def __str__(self):
        """ print out defined attributes of a CSS_Event object """
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

    def parse(self, line):
        """
        Parse a line of a CSS Event table and populate a CSS_Event object
        """

        try:
            # generate fields for this event by parsing line
            ptr = 0
            ptr, self.evid = parse_attribute(line, 8, ptr, 'i') 
            ptr, self.evname = parse_attribute(line, 15, ptr+1, 's') 
            ptr, self.evid = parse_attribute(line, 8, ptr, 'i') 
            ptr, self.evname = parse_attribute(line, 15, ptr+1, 's') 
            self.prefor = int(line[23:31])
            self.auth = line[31:46].strip().decode()
            self.commid = int(line[46:54])
            self.lddate = float(line[54:])

            # check for nulls
            if self.evid == NULLID:
                self.evid = None
            if self.evname == NULLSTRING:
                self.evname = None
            if self.prefor == NULLID:
                self.prefor = None
            if self.auth == NULLSTRING:
                self.auth = None
            if self.commid == NULLID:
                self.commid = None
            if self.lddate == NULLTIME:
                self.lddate = None
        except:
            print("Error reading event from line: %s" % line)

class CSS_Origin(object):
    """
    Blueprint for a CSS_Origin object

    The format of an Antelope/CSS3.0 schema origin table is taken from $ANTELOPE/data/schemas/css3.0
    (part of the commercial Antelope system)

    The fields are:

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

    def __init__(self, lat=None, lon=None, depth=None, time=None, orid=1, evid=None, 
        jdate=None, ndef=None, nass=None, 
        ndp=None, grn=None, srn=None, etype=None, review=None, mbid=None, mb=None, 
        msid=None, ms=None, mlid=None, ml=None, 
        depdp=None, dtype=None, algorithm=None, auth=None, commid=None, lddate=None):
        """ initialise an empty Origin object """
        if VERBOSE:
            print("Creating CSS Origin object")
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
        """ print out defined attributes of a CSS_Origin object """
        str = "CSS Origin object with attributes:\n"
        if self.lat:
            str += "\tlat: %s\n" % self.lat
        if self.lon:
            str += "\tlon: %s\n" % self.lon
        if self.depth:
            str += "\tdepth: %s\n" % self.depth
        if self.time:
            str += "\ttime: " + self.time.__str__() + "\n"
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
            str += "\tlddate: " + self.lddate.__str__() + "\n"
        if self.arrivals:
            str += "\tarrivals: %f\n" % len(self.arrivals)
        if self.netmags:
            str += "\tnetmags: %f\n" % len(self.netmags)
        return str

#    def add_arrival(self, arrivalobject):
#        self.arrivals.append(arrivalobject)

#    def add_netmag(self, netmagobject):
#        self.netmags.append(netmagobject)

#    def add_origerr(self, origerrobject):
#        self.origerr = origerrobject 

    def parse(self, line):
        """
        Parse a line of a CSS Origin table and return a CSS_Origin object
        """

        try:
            # generate fields for this origin by parsing line
            ptr = 0
            ptr, self.lat = parse_attribute(line, 9, ptr, 'f') 
            ptr, self.lon = parse_attribute(line, 9, ptr+1, 'f') 
            ptr, self.depth = parse_attribute(line, 9, ptr+1, 'f') 
            ptr, self.time = parse_attribute(line, 17, ptr+1, 't') 
            ptr, self.orid = parse_attribute(line, 8, ptr+1, 'i')
            ptr, self.evid = parse_attribute(line, 8, ptr+1, 'i')
            ptr, self.jdate = parse_attribute(line, 8, ptr+1, 'i')
            ptr, self.nass = parse_attribute(line, 4, ptr+1, 'i')
            ptr, self.ndef = parse_attribute(line, 4, ptr+1, 'i')
            ptr, self.ndp = parse_attribute(line, 4, ptr+1, 'i')
            ptr, self.grn = parse_attribute(line, 8, ptr+1, 'i')
            ptr, self.srn = parse_attribute(line, 8, ptr+1, 'i')
            ptr, self.etype = parse_attribute(line, 2, ptr+1, 's')
            ptr, self.review = parse_attribute(line, 4, ptr+1, 's')
            ptr, self.depdp = parse_attribute(line, 9, ptr+1, 'f')
            ptr, self.dtype = parse_attribute(line, 1, ptr+1, 's')
            ptr, self.mb = parse_attribute(line, 7, ptr+1, 'f')
            ptr, self.mbid = parse_attribute(line, 8, ptr+1, 'i')
            ptr, self.ms = parse_attribute(line, 7, ptr+1, 'f')
            ptr, self.msid = parse_attribute(line, 8, ptr+1, 'i')
            ptr, self.ml = parse_attribute(line, 7, ptr+1, 'f')
            ptr, self.mlid = parse_attribute(line, 8, ptr+1, 'i')
            ptr, self.algorithm = parse_attribute(line, 15, ptr+1, 's')
            ptr, self.auth = parse_attribute(line, 15, ptr+1, 's')
            ptr, self.commid = parse_attribute(line, 8, ptr+1, 'i')
            ptr, self.lddate = parse_attribute(line, 17, ptr+1, 't')
        except:
            print("Error parsing origin from line\n: %s" % line)
        return

if __name__ == "__main__":

    # read in the Antelope demo database into a CSS_Catalog object
    c = CSS_Catalog()
    c.read('/opt/antelope/data/db/demo/demo')

    # Inspect it
    print(c.events[0])
    print(c.events[0].origins[0])

    # export the CSS_Catalog object to an ObsPy Catalog object
    c2 = c.to_obspy()
    print(c2)
    print(c2.events[0])
    print(c2.events[0].origins[0])

    # write the ObsPy Catalog object to QuakeML
    c2.write("/tmp/my_custom_events.xml", format="QUAKEML")

