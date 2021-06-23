#! /usr/bin/env python
""" arrival class, Glenn Thompson 2010/05/01 """
    
    
class Arrival(object):
    
    def __init__(self):
        """ initialise an empty Arrival object """
    self= numpy.array([])

#    def __str__(self):
#        """ print out the attributes of an event object """
#   print self
 
    def import_random(self, n):
        """ create N random events """
        dtnow = datetime.datetime.utcnow()
    arid = numpy.arange(n)
    etime = converttime.utnow() - numpy.random.random(n) * 86400
    sta = numpy.array(['DUMM']*n)
    chan = numpy.array(['EHZ']*n])
    net = numpy.array(['AV']*n])
    phase = numpy.array(['P']*n])
    self = numpy.rec.fromarrays([arid, etime, sta, chan, net, phase],  names = 'arid, etime, sta, chan, net, phase' )

    def import_antelope(self, estart, eend, dbname, region, expr):
        """ load events from an Antelope database """
        """ not tried yet, and does not support day/month archives """
    evid = etype = orid = etime = lon = lat = depth = ml = ms = mb = []
        db = dbopen(dbname, 'r')
        dbe = dbopen_table(db, 'event')
        dbo = dbopen_table(db, 'origin')
    dbj = dbjoin(dbe, dbo)
    dbs = dbsubset(dbj, 'orid == prefor')   
    dbas = dbopen_table(db, 'assoc')
    dbj = dbjoin(dbs, dbas)
    dbar = dbopen_table(db, 'arrival')  
    dbj = dbjoin(dbj, dbar)
        dbs = dbsubset(dbs, expr)
        nrecs = dbquery(dbs, 'dbRECORD_COUNT')
    # Try vectorized version ?
    [arid, etime, sta, chan, net, phase] = dbgetv(dbd, 'arid', 'etime', 'sta', 'chan', 'net', 'phase')
        dbclose(db)
    self = numpy.rec.fromarrays([arid, etime, sta, chan, net, phase],  names = 'arid, etime, sta, chan, net, phase' )
    
        
    def import_antelope_wrapper(self, filepattern, startdate, enddate, expr):
        """ wrapper for loadaef over multiple months """
        # loop over day/month archives
            # generate db name
            # import_antelope   

if __name__ == "__main__":
    a = Arrival()
    n = input('Number of arrivals ? ')
    a.import_random(n)
    print a.__str__()    
    

