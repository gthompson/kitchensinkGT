"""
A Python reader for Antelope CSS event databases without using Antelope/Datascope 
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
import struct
import numpy as np
from obspy import UTCDateTime, Trace, Stream

def readCSS(dbname, **kwargs):
    """
    Reads an Antelope CSS3.0 event database and returns a catalog object.

    :type dbname: string
    :param dbname: CSS3.0 database to be read.
    :rtype: :class:`~obspy.core.event.Catalog`
    :returns: Catalog with Events specified by given database name.

    A CSS3.0 event database consists of the following tables:
    - event
    - origin (1 or many per event)
    - netmag (0, 1 or many per origin)
    - origerr (0 or 1 per origin)
    - assoc (links origin and arrival)
    - arrival (each arrival may be associated with 0, 1 or many origins)

    """

    """
    READ IN THE EVENT TABLE - this is what it looks like in Antelope/CSS3.0 schema:
    Relation event
        Fields ( evid evname prefor auth commid lddate )
        evid [0-7] - integer(%8ld), null -1, evid>0
        evname [8-22] - string(%-15s), null -
        prefor [23-30]  - integer(%8ld), null -1, prefor>0
        auth [31-45] - string(%-15s), null -
        commid [46-53] - integer(%8ld), null -1, commid>0
        lddate [54-70] - time(%17.5f), null -9999999999.99900, epoch seconds
    """
    # initialize
    c = CSSCatalog()
    
    # check if event table exists
    if(os.path.exists(os.path.join(dbname, '.event'))):
        with open(filename, "r") as fh:
            lines = fh.readlines()

    # process line by line
    line_num = 0
        for line in lines:
        try:
            line_num += 1
            e = CSSEvent()
            e.evid = int(line[0:8])
            e.evname = line[8:23].strip().decode()
            e.prefor = int(line[23:31])
            e.auth = line[31:46].strip().decode()
            e.commid = int(line[46:54])
            e.lddate = float(line[54:])
            c.events.append(e)
        except:
            print("Error reading event at line %d" % line_num)

    """
    READ IN THE ORIGIN TABLE - this is what it looks like in Antelope/CSS3.0 schema:
    Relation origin
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

    # check if origin table exists
    if(os.path.exists(os.path.join(dbname, '.origin'))):
        with open(filename, "r") as fh:
            lines = fh.readlines()

    # process line by line
    line_num = 0
        for line in lines:
        try:
            line_num += 1
            o = CSSOrigin()
            o.lat = float(line[0:9])
            o.lon = float(line[9:18])
            o.depth = float(line[18:27])
            o.time = float(line[27:44])
            o.orid = int(line[44:52])
            o.evid = int(line[52:60])
            o.jdate = int(line[60:68])
            o.nass = int(line[68:72])
            o.ndef = int(line[72:76])
            o.ndp = int(line[76:80])
            o.grn = int(line[80:88])
            o.srn = int(line[88:96])
            o.etype = line[96:98].strip().decode()
            o.review = line[98:102].strip().decode()
            o.depdp = float(line[102:111])
            o.dtype = line[102:103].strip().decode()
            o.mb = float(line[103:110])
            o.mbid = int(line[110:118])
            o.ms = float(line[118:125])
            o.msid = int(line[125:133])
            o.ml = float(line[133:140])
            o.mlid = int(line[140:148])
            o.algorithm = line[148:163].strip().decode()
            o.auth = line[163:178].strip().decode()
            o.commid = int(line[178:186])
            o.lddate = float(line[186:203])
            c.origins.append(o)

        except:
            print("Error reading origin at line %d" % line_num)

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

    return c
