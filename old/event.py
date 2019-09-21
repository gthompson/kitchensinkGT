#! /usr/bin/env python
""" event class, Glenn Thompson 2010/05/01 """
    
class event(object):

    def __init__(self, evid=-1, etype='u'):
        """ initialise an empty Event object """
        self.evid = evid
        self.etype = etype
        self.origin = {}
        self.arrival = []

    def __str__(self):
        """ print out the attributes of an event object """
        print '%5d' % self.evid,', %s' % self.etype, ', %s' % self.origin, ', %s' % self.arrival

    def add_origin(self, originobject):
        self.origin = originobject

    def add_arrival(self, arrivalobject):
        self.arrival.append(arrivalobject)

if __name__ == "__main__":
    import sys
    sys.path.append('F:/avoseis')
    import origin
    import arrival
    #import location
    e = event(1, 'r')
    #o = origin.origin()
    #a = arrival.arrival()
    #e.add_origin(o)
    #e.add_arrival(a)
    print e

