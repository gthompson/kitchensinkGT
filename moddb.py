import antelope.datascope as datascope
import matplotlib as mpl
import numpy as np
import glob

def remove_database(dbpath, ask=True):
    if os.path.exists("%s" % dbpath):
        doit = True
        if ask:
            doit = yn_choice("%s database exists. Remove?" % dbpath)
        if doit:
            os.remove(dbpath)
            for dbfile in glob.glob("%s.*" % dbpath):
                os.remove(dbfile)
            print "database %s removed" % dbpath
        else:
            print "You have chosen not to remove %s" % dbpath
    else:
        print "database %s not found" % dbpath

def dbsubset2db(dbpath, subset_expr, subsetdbpath):
    # open the origin table, join to event table, subset for preferred origins
    db = datascope.dbopen( dbpath, 'r')
    dborigin = db.lookup( table = 'origin' )
    dborigin = dborigin.join('event')
    dbnetmag = db.lookup( table = 'netmag' )
    dborigin = dborigin.subset("orid == prefor")
    print "subset %s with %s" % (dbpath, subset_expr)
    dborigin = dborigin.subset(subset_expr)
    dborigin = dborigin.sort('time')
    n = dborigin.nrecs()
    if n>0:    
        # ask to remove database if it already exists
        remove_database(subsetdbpath, False)
        print "unjoin to %s" % subsetdbpath
        dborigin.unjoin(subsetdbpath)
    else:
        print "no records to save to new database\n"
    db.free()

def dbgetorigins(dbpath, subset_expr):
    # open the origin table, join to event table, subset for preferred origins
    db = datascope.dbopen( dbpath, 'r')
    dborigin = db.lookup( table = 'origin' )
    dborigin = dborigin.join('event')
    dbnetmag = db.lookup( table = 'netmag' )
    dborigin = dborigin.subset("orid == prefor")

    # apply the optional subset expression if there is one, order by time, and display number of events.
    if subset_expr:
        dborigin = dborigin.subset(subset_expr)
    dborigin = dborigin.sort('time')
    n = dborigin.nrecs()
    #print "- number of events = {}".format(n)

    # if size of arrays already known, preallocation much faster than recreating each time with append
    dictorigin = dict()
    origin_id = np.empty(n)
    origin_ml = np.empty(n)
    origin_epoch = np.empty(n)

    # load origins from database and store them in a dictionary
    for dborigin[3] in range(n):
        (origin_id[dborigin[3]], origin_ml[dborigin[3]], origin_epoch[dborigin[3]]) = dborigin.getv('orid','ml','time')
        if origin_ml[dborigin[3]] < -1.0:
            db2 = dbnetmag.subset("orid == %d" % origin_id[dborigin[3]])
            maxmag = -1.0
            n_netmag = db2.nrecs()
            if n_netmag > 0:
                for db2[3] in range(n_netmag):
                    (magtype, magnitude) = db2.getv('magtype', 'magnitude')
                    if magnitude>maxmag:
                        maxmag = magnitude
            origin_ml[dborigin[3]] = maxmag
            
    dictorigin['id'] = origin_id
    dictorigin['ml'] = origin_ml
    dictorigin['time'] = mpl.dates.epoch2num(origin_epoch)

    # close the database and free the memory. 
    # It seems that db.close and db.free both close the database, and closing twice produces error
    db.free()

    return dictorigin, n

def readplacesdb(dbplacespath):
    """
    should be able to do this with code like:
    
    db = datascope.dbopen( dbpath, 'r')
    db = db.lookup( table = 'places' )
    db = db.sort('lon')
    n = db.nrecs()
    if n > 0:
        print "- number of places = {}".format(n)
        for db[3] in range(n):
                (placename, placetype, placelat, placelon, placeelev) = db.getv('place','placetype','lat','lon','elev')
    
     but this does not work. Neither does the Perl equivalent. So although dbe opens the places database, there might
     be a problem with the database. So just treat the table as a file and parse it instead.
    """
    dbtablepath = dbplacespath + ".places"
    placename = list()
    placetype = list()
    placelat = list()
    placelon = list()
    placeelev = list()
    placeradius = list()
    if os.path.exists(dbtablepath):
        f = open(dbtablepath)
        lines = f.readlines()
        f.close()
        for line in lines:
            elements = line.split()
            placename.append(elements[4])
            placetype.append(elements[5])
            placelat.append(elements[0])
            placelon.append(elements[1])
            placeelev.append(elements[2])
            placeradius.append(elements[3])
    else:
        print dbtablepath + " does not exist"        
    return {'place':placename, 'placetype':placetype, 'lat':placelat, 'lon':placelon, 'elev':placeelev, 'radius':placeradius}

def read_volcanoes(): # THIS IS A PREVIOUS VERSION I WROTE FOR GETTING VOLCANOES FROM PLACESDB THAT SEEMS TO WORK
    dbplacespath = 'places/volcanoes'
    dbhandle = datascope.dbopen( dbplacespath, 'r')
    dbptr = dbhandle.lookup( table = 'places' )
    n = dbptr.nrecs()
    dictplaces = dict()
    for dbptr[3] in range(n):
        thisrecord = {'place': "%s" % dbptr.getv('place'), 'lat': "%s" % dbptr.getv('lat'), 'lon': "%s" % dbptr.getv('lon') }
        dictplaces[dbptr[3]] =  thisrecord
    dbhandle.free()
    dbhandle.close()
    return dictplaces

