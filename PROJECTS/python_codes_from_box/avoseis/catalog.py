#! /usr/bin/env python
""" Catalog class, Glenn Thompson 2010/04/29 """
""" This is intended to be a full representation of the CSS3.0
    event, origin, arrival and assoc schema (prefors only)
    The primary purpose is to provide conversion routines between
    different catalog types """
    
class Catalog(object):
    
    def __init__(self):
        """ initialise an empty Catalog object """
        """ a Catalog object consists of:
            event.evid
            event.etype
            event.origin.orid
            event.origin.datetime
            event.origin.location.lon
            event.origin.location.lat
            event.origin.location.depth
            event.origin.ml
            event.origin.mb
            event.origin.ms
            event.origin.mw
            event.arrival[j].sta
            event.arrival[j].chan
            event.arrival[j].net
            event.arrival[j].phase
            event.arrival[j].datetime

            So the object types are:
                event
                origin
                arrival
                location
                datetime
            
        """
            
        self = []

    def __str__(self):
        """ print out the attributes of an event object """
        print 'num\tevid\tlon\tlat\tdepth\tmag'
        for i in range(len(self.evid)):
            #print i,'\t',self.evid[i],'\t',self.lon[i],'\t',self.lat[i],'\t: ',self.depth[i],'\t',self.mag[i]
            print '%4d' % i,' %5d' % self.evid[i],' %7.2f' % self.lon[i],' %7.2f' % self.lat[i],' %5.1f' % self.depth[i],' %4.1f' % self.mag[i]

    def add_event(self, event):
        self.append(event)
            
    def import_random(self, n):
        """ create N random events """
        import random
        import datetime
        import converttime
        d = datetime.datetime.utcnow()
        for i in range(n):
            epochtime = convertime.datetime2epoch(d.utcnow() - datetime.timedelta(random.random(),0,0))
            lon = (random.random() * 360 - 180)
            lat = (random.random() * 180 - 90)
            depth = (random.random() * 30)
            ml = (random.random() * 3 - 0.5)
            etype = 'l'
            o = origin.origin(i, epochtime, lon, lat, depth, ml)
            e = event(i, etype)
            e.add_origin(o)

    def import_antelope(self, file, expr):
        """ load events from an Antelope database """
        """ not tried yet, and does not support day/month archives """
        import converttime
        db = dbopen(file, 'r')
        dbe = dbopen_table(db, 'event')
        dbo = dbopen_table(db, 'origin')
        dbj = dbjoin(dbe, dbo)
        dbs = dbsubset(dbj, 'orid == prefor')
        dbs = dbsubset(dbs, expr)
        nrecs = dbquery(dbs, 'dbRECORD_COUNT')
        for i in range(nrecs):
            dbs[3] = i-1
            [etime, evid, etype, orid, lon, lat, depth, ml, ms, mb] = dbgetv(dbs, 'time', 'evid', 'etype', 'orid', 'lon', 'lat', 'depth', 'ml', 'ms', 'mb')
            o = origin.origin(orid, etime, lon, lat, depth, ml, ms, mb)
            e[i] = event(evid, etype)
            e[i].add_origin(o)
        dbclose(db)
        
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

    def import_hypo71(self, file)
        """ load event data from a hypo71 catalog """

    def import_hypoellipse(self, file)
        """ load event data from a hypoellipse catalog """

    def import_hypoinverse(self, file)
        """ load event data from a hypoinverse catalog """

    def import_hypocenter(self, file)
        """ load event data from a hypocenter catalog """

    def plotlonlat(self, FontSize):
        """ plot latitude versus longitude (Map view) """
        import pylab as plt
        plt.plot(self.lon, self.lat, 'o')
        #plt.title('Map view')
        plt.xlabel('Longitude', fontsize=FontSize)
        plt.ylabel('Latitude', fontsize=FontSize)
        #ah = plt.gca()
        #xtl = ah.get_xticklabels()
        #print xtl
        #ah.set_xticklabels(xtl, fontsize=FontSize)

    def plotdepthtime(self, FontSize):
        """ plot depth against time  """
        import pylab as plt
        plt.plot_date(self.time, self.depth,'o')
        #plt.title('Depth-time view')
        plt.xlabel('time', fontsize=FontSize)
        plt.ylabel('depth', fontsize=FontSize)
        plt.gca().invert_yaxis()

    def plotmagtime(self, FontSize):
        """ plot magnitude versus time """
        import pylab as plt
        plt.plot_date(self.time, self.mag)
        #plt.title('Mag-time view')
        plt.xlabel('time', fontsize=FontSize)
        plt.ylabel('magnitude', fontsize=FontSize)

    def plotdepthlon(self, FontSize):
        """ plot depth vs. longitude """
        import pylab as plt
        plt.plot(self.lon, self.depth,'o')
        #plt.title('Depth-lon view')
        plt.xlabel('longitude', fontsize=FontSize)
        plt.ylabel('depth', fontsize=FontSize)
        plt.gca().invert_yaxis()

    def plotlatdepth(self, FontSize):
        """ plot latitude vs. depth """
        import pylab as plt
        plt.plot(self.depth, self.lat,'o')
        #plt.title('Lat-depth view')
        plt.xlabel('depth', fontsize=FontSize)
        plt.ylabel('latitude', fontsize=FontSize)

    def volplot(self, pngfile):
        """ Replicate Guy's VolPlot routines """
        import pylab as plt
        FontSize = 12
        plt.close('all')
        plt.figure(1)
        plt.axes([0.1, 0.45, 0.5, 0.5])
        self.plotlonlat(FontSize)
        plt.axes([0.1, 0.25, 0.5, 0.15])
        self.plotdepthlon(FontSize)
        plt.axes([0.1, 0.05, 0.8, 0.15])    
        self.plotdepthtime(FontSize)
        plt.axes([0.7, 0.45, 0.2, 0.5])    
        self.plotlatdepth(FontSize)
        #plt.show()
        plt.savefig(pngfile)

    def datebin(self, binsize):
        """ bin the data based on a bin size of binsize days """
        import converttime
        #for i in range(len(self.evid)):
        #    etime[i] = converttime.epoch2datetime(self.time[i])
        etime = converttime.epoch2datetime(self.time)

    def add_epochtime(self):
        """ add epoch time to the event object """
        import converttime
        #for i in range(len(self.evid)):
        #    etime[i] = converttime.datetime2epoch(self.time[i])
        self.etime = converttime.datetime2epoch(self.time)

    def export_xml(self):
        """ print out event object as xml """
        for i in range(len(self.evid)):
            print '<event>'
            print '\t<evid>',self.evid[i],'</evid>'
            print '\t<time>',self.time[i],'</time>' # convert datetime to ?      
            print '\t<lon>',self.lon[i],'</lon>'
            print '\t<lat>',self.lat[i],'</lat>'
            print '\t<depth>',self.depth[i],'</depth>'
            print '\t<mag>',self.mag[i],'</mag>'
            print '</event>'

    def export_kml(self):
        """ print out event object as kml """
        for i in range(len(self.evid)):
            print '<event>'
            print '\t<evid>',self.evid[i],'</evid>'
            print '\t<time>',self.time[i],'</time>' # convert datetime to GoogleEarth time?        
            print '\t<lon>',self.lon[i],'</lon>'
            print '\t<lat>',self.lat[i],'</lat>'
            print '\t<depth>',self.depth[i],'</depth>'
            print '\t<mag>',self.mag[i],'</mag>'
            print '</event>'

    def export_antelope(self, dbname):
        """ print out event object as an antelope origin table """
        db = dbopen(dbname, 'a+')
        dbe = dbopen_table(db + '.event', 'a+')
        dbo = dbopen_table(db + '.origin', 'a+')
        for i in range(len(self.evid)):
            # create a new row?
            # use dbaddv or dbputv ?
            # origin table
            dbputv(dbo, 'orid', self.evid[i])
            dbputv(dbo, 'time', self.time[i]) # convert datetime to epoch time
            dbputv(dbo, 'lon', self.lon[i])
            dbputv(dbo, 'lat', self.lat[i])
            dbputv(dbo, 'depth', self.depth[i])
            dbputv(dbo, 'ml', self.mag[i])
            # event table
            dbputv(dbe, 'evid', self.evid[i])
            dbputv(dbe, 'prefor', self.evid[i])
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
    print e.__str__()    
    #volplot('F:/volplot.png')


    


    
