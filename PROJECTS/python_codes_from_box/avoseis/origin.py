#! /usr/bin/env python
""" origin class, Glenn Thompson 2010/04/29 """

class origin(object):
    
    import location
    
    def __init__(self, orid=-1, epochtime=0, lon=0, lat=0, depth=0, ml=-999.9, mb=-999.9, ms=-999.9, mw=-999.9):
        """ initialise an origin object """
        import converttime
        self.orid = orid
        self.datetime = converttime.epoch2datetime(epochtime)
        self.location = location.location(lon, lat, depth, 'down')
        self.ml = ml
        self.mb = mb
        self.ms = ms
        self.mw = mw

    def __str__(self):
        """ print out the attributes of an origin object """
        #print 'orid = %d' % self.orid, ', time = %s' % self.datetime, ', %s' % self.location, ', ml = %4.1f' % self.ml
        print 'orid = %d' % self.orid, ', time = %s' % self.datetime, ', ml = %4.1f' % self.ml           

if __name__ == "__main__":
    import random
    l = [origin(), origin(2, 1200000000, random.random() * 360 - 180, random.random() * 180 - 90, random.random() * 30, 2.5)]
    for i in range(len(l)):
        print l[i].__str__()
    print l[1]
    

