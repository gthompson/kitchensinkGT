default_float = -9999.9999
default_integer = -1
default_string = '-'

class arrival: # from tables arrival, assoc, stamag
    def __init__(self, atime, iphase):
        self.atime = atime
        self.iphase = iphase
        self.stamag = stamag

class event: # from tables event, origin, netmag
    def __init__(self, time=default_float, lon=default_float, lat=default_float, depth=default_float, mag=default_float, magtype=default_string, nass=default_integer, ndef=default_integer, etype=default_string, auth=default_string, arrivals=[] ):
        self.time = time
        self.lon = lon 
        self.lat = lat
        self.depth = depth
        self.mag = mag
        self.magtype = magtype
        self.nass = nass
        self.nass = ndef
        self.etype = etype
        self.auth = auth
        self.arrivals = arrivals # call eventobj.loadarrivals() to populate

class catalog: # collection of events
    def __init__(self, dbpath, subset_expr=default_string, get_arrivals=False):
        self.dbpath = dbpath
        self.subset_expr = subset_expr
        self.events = []
        # load events from database
        print("Load events from %s " % self.dbpath)
        if get_arrivals:
            # load arrivals too
            print("Load arrivals too")


    
    
