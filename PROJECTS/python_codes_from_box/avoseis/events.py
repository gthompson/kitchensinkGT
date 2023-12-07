#! /usr/bin/env python
""" Event class, Glenn Thompson 2010/04/29 """

# might eventually just want to import specific modules used to shorten load time
import numpy
import random
import datetime
import converttime
import scipy
import pylab
vector_epoch2datetime = scipy.vectorize(converttime.epoch2datetime)

class Event(object):
    
    def __init__(self):
        """ initialise an empty Event object """
	self = numpy.array([])

    #def __str__(self):
    #    """ print out the attributes of an event object """
	#print self
        ##print 'num evid time                                        lon        lat       depth  mag'
        ##for i in range(len(self)):
        ##    #print i,'\t',self.evid[i],'\t',self.lon[i],'\t',self.lat[i],'\t: ',self.depth[i],'\t',self.mag[i]
        ##    print '%4d' % i,' %5d' % self.evid[i],' %s' % self.time[i],' %7.2f' % self.lon[i],' %7.2f' % self.lat[i],' %4.1f' % self.depth[i],' %3.1f' % self.mag[i]
            
    def import_random(self, n):
        """ create N random events """
        dtnow = datetime.datetime.utcnow()
	evid = etype = orid = etime = lon = lat = depth = ml = ms = mb = []
	evid = orid = numpy.arange(n)
	etype = numpy.array(['r']*n)
	orid = numpy.arange(n)
	etime = converttime.utnow() - numpy.random.random(n) * 86400
	lon = numpy.random.random(n) * 360 - 180
	lat = numpy.random.random(n) * 180 - 90
	depth = numpy.random.random(n) * 30
	ml = numpy.random.random(n) * 3 - 0.5
	ms = mb = numpy.ones(n) * -999.9
	self = numpy.rec.fromarrays([evid, etype, orid, etime, lon, lat, depth, ml, ms, mb],  names = 'evid, etype, orid, etime, lon, lat, depth, ml, ms, mb' )

    def import_antelope(self, estart, eend, dbname, region, expr):
        """ load events from an Antelope database """
        """ not tried yet, and does not support day/month archives """
	evid = etype = orid = etime = lon = lat = depth = ml = ms = mb = []
        db = dbopen(dbname, 'r')
        dbe = dbopen_table(db, 'event')
        dbo = dbopen_table(db, 'origin')
        dbj = dbjoin(dbe, dbo)
        dbs = dbsubset(dbj, 'orid == prefor')
        dbs = dbsubset(dbs, expr)
        nrecs = dbquery(dbs, 'dbRECORD_COUNT')
	# Try vectorized version ?
	[evid, etype, orid, etime, on, lat, depth, ml, mb, ms] = dbgetv(dbs, 'evid', 'etype', 'orid', 'time', 'lon', 'lat', 'depth', 'ml', 'mb', 'ms')
        # Otherwise use this
	#for i in range(nrecs):
        #    dbs[3] = i-1
        #    [evid, etype, orid, etime, on, lat, depth, ml, mb, ms] = dbgetv(dbs, 'evid', 'etype', 'orid', 'time', 'lon', 'lat', 'depth', 'ml', 'mb', 'ms')
        #    self.etime.append(etime)
        #    self.datetime.append(converttime.epoch2datetime(etime)) # convert from epoch time to datetime
        #    self.evid.append(evid)
        #    self.lon.append(lon)
        #    self.lat.append(lat)
        #    self.depth.append(depth)
	#    self.location = location.location(lon, lat, depth, 'down')
        #    self.mag.append(mag)
        dbclose(db)
	self = numpy.rec.fromarrays([evid, etype, orid, etime, lon, lat, depth, ml, ms, mb], names = 'evid, etype, orid, etime, lon, lat, depth, ml, ms, mb' )
	
        
    def import_antelope_wrapper(self, filepattern, startdate, enddate, expr):
        """ wrapper for loadaef over multiple months """
        # loop over day/month archives
            # generate db name
            # import_antelope

    def import_aef(self, file, expr):
        """ load events from an MVO aef month summary """
        # open file
        # loop
            # read line
            # parse line
            # convert date/time to datetime
            # append items to lists
        # close file
        
    def import_aef_wrapper(self, filepattern, startdate, enddate):
        """ wrapper for loadaef over multiple months """
        # loop over months
            # generate file name
            # import_aef

    def import_seisan(self, filepattern, startdate, enddate):
        """ load event data from a Seisan REA database """
        # loop over year/month directories
            # read file names
            # loop over files
                # parse event/origin/arrival information
                # append to event lists (time, lat, lon, ...)

    def import_hypo71(self, file):
        """ load event data from a hypo71 catalog """

    def import_hypoellipse(self, file):
        """ load event data from a hypoellipse catalog """

    def import_hypoinverse(self, file):
        """ load event data from a hypoinverse catalog """

    def import_hypocenter(self, file):
        """ load event data from a hypocenter catalog """
		
    def plotlonlat(self, FontSize):
        """ plot latitude versus longitude (Map view) """
        pylab.plot(self['lon'], self['lat'], 'o')
        #pylab.title('Map view')
        pylab.xlabel('Longitude', fontsize=FontSize)
        pylab.ylabel('Latitude', fontsize=FontSize)
        #ah = pylab.gca()
        #xtl = ah.get_xticklabels()
        #print xtl
        #ah.set_xticklabels(xtl, fontsize=FontSize)

    def plotdepthtime(self, FontSize):
        """ plot depth against time  """
        pylab.plot_date(vector_epoch2datetime(self['etime']), self['depth'],'o')
        #pylab.title('Depth-time view')
        pylab.xlabel('time', fontsize=FontSize)
        pylab.ylabel('depth', fontsize=FontSize)
        pylab.gca().invert_yaxis()

    def plotmagtime(self, FontSize):
        """ plot magnitude versus time """
        pylab.plot_date(vector_epoch2datetime(self['etime']), self['ml'])
        #pylab.title('Mag-time view')
        pylab.xlabel('time', fontsize=FontSize)
        pylab.ylabel('magnitude', fontsize=FontSize)

    def plotdepthlon(self, FontSize):
        """ plot depth vs. longitude """
        pylab.plot(self['lon'], self['depth'],'o')
        #pylab.title('Depth-lon view')
        pylab.xlabel('longitude', fontsize=FontSize)
        pylab.ylabel('depth', fontsize=FontSize)
        pylab.gca().invert_yaxis()

    def plotlatdepth(self, FontSize):
        """ plot latitude vs. depth """
        pylab.plot(self['depth'], self['lat'],'o')
        #pylab.title('Lat-depth view')
        pylab.xlabel('depth', fontsize=FontSize)
        pylab.ylabel('latitude', fontsize=FontSize)

    def volplot(self, pngfile):
        """ Replicate Guy's VolPlot routines """
        FontSize = 12
        pylab.close('all')
        pylab.figure(1)
	print 'Plotting lat vs.lon'
        pylab.axes([0.1, 0.45, 0.5, 0.5])
        self.plotlonlat(FontSize)
	print 'Plotting lon vs. depth'
        pylab.axes([0.1, 0.25, 0.5, 0.15])
        self.plotdepthlon(FontSize)
	print 'Plotting depth vs. time'
        pylab.axes([0.1, 0.05, 0.8, 0.15])    
        self.plotdepthtime(FontSize)
	print 'Plotting lat vs. depth'
        pylab.axes([0.7, 0.45, 0.2, 0.5])    
        self.plotlatdepth(FontSize)
        #pylab.show()
	print 'Saving figure'
        pylab.savefig(pngfile)

    def plotcounts(self, binsize, pngfile):
        """ bin the data based on a bin size of binsize days """
        #vector_datetime2epoch = scipy.vectorize(converttime.datetime2epoch)
        #etime = vector_d2e(self['etime'])
	etime = self['etime']
	etimemin = scipy.floor(min(etime)/86400) * 86400
	etimemax = scipy.ceil(max(etime)/86400) * 86400
	r = [etimemin, etimemax]
	nbins = (etimemax - etimemin) / (binsize * 86400)
	h = scipy.histogram(etime, nbins, r)
	ecount = h.__getitem__(0)
	ebin = h.__getitem__(1)
	#vector_e2d = scipy.vectorize(converttime.epoch2datetime)
	edt = vector_epoch2datetime(ebin)
	pylab.plot_date(edt, ecount, linestyle='steps-mid')
	pylab.savefig(pngfile)

    def add_epochtime(self):
        """ add epoch time to the event object """
        self.etime = converttime.datetime2epoch(self.time)

    def export_xml(self):
        """ print out event object as xml """
        #for i in range(len(self)):
        #    print '<event>'
        #    print '\t<evid>',self.evid[i],'</evid>'
        #    print '\t<time>',self.time[i],'</time>' # convert datetime to ?      
        #    print '\t<lon>',self.lon[i],'</lon>'
        #    print '\t<lat>',self.lat[i],'</lat>'
        #    print '\t<depth>',self.depth[i],'</depth>'
        #    print '\t<mag>',self.mag[i],'</mag>'
        #    print '</event>'

    def export_kml(self):
        """ print out event object as kml """
        #for i in range(len(self.evid)):
        #    print '<event>'
        #    print '\t<evid>',self.evid[i],'</evid>'
        #    print '\t<time>',self.time[i],'</time>' # convert datetime to GoogleEarth time?        
        #    print '\t<lon>',self.lon[i],'</lon>'
        #    print '\t<lat>',self.lat[i],'</lat>'
        #    print '\t<depth>',self.depth[i],'</depth>'
        #    print '\t<mag>',self.mag[i],'</mag>'
        #    print '</event>'

    def export_antelope(self, dbname):
        """ print out event object as an antelope origin table """
        db = dbopen(dbname, 'a+')
        dbe = dbopen_table(db + '.event', 'a+')
        dbo = dbopen_table(db + '.origin', 'a+')
        #for i in range(len(self)):
            # create a new row?
            # use dbaddv or dbputv ?
            # origin table
        #    dbputv(dbo, 'orid', self.evid[i])
        #    dbputv(dbo, 'time', self.time[i]) # convert datetime to epoch time
        #    dbputv(dbo, 'lon', self.lon[i])
        #    dbputv(dbo, 'lat', self.lat[i])
        #    dbputv(dbo, 'depth', self.depth[i])
        #    dbputv(dbo, 'ml', self.mag[i])
        #    # event table
        #    dbputv(dbe, 'evid', self.evid[i])
        #    dbputv(dbe, 'prefor', self.evid[i])
        dbclose(db)

    def export_seisan(self, filepattern):
        """ export event object as Seisan REA files in a year/month structure """

    def event_gui(self):
        """ tkinter gui for interactive importing, plotting and exporting of event objects """

    def eev(self):
        """ catalog tool based on Seisan's eev """
        """ or it could be a tk interface with a window showing info about this event, plus prev and next arrows """
        """ plus a tool to leap based on yyyymmdd[hh[mm]] """

if __name__ == "__main__":
    e = Event()
    n = input('Number of events ? ')
    e.import_random(n)
    #print e.__str__()    
    
    #e.plotcounts(1.0/24, 'F:/hist.png')
    e.volplot('F:/volplot.png')
    


    
