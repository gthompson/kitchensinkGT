AVOSEIS=os.getenv("AVOSEIS");
#sys.path.append('~/src/python_gt/AVOSEIS_PYTHON')
sys.path.append(AVOSEIS + "/bin")
import antelope.datascope as datascope
import matplotlib as mpl
if 'DISPLAY' in os.environ.keys():
        mpl.use("Agg")
#import numpy as np
import matplotlib.pyplot as plt
import modgiseis
#import time
import datetime
import getopt
import dbploteventrate
######################################################

def usage():
        print 'Usage: '+sys.argv[0]+' <catalogpath> <dbplacespath> <outputdir>'
        print """
        <catalogpath> must have an origin and event table present
        <dbplacespath> is a list of volcanoes and their lat, lon, elev and radii in places_avo_1.3 schema format
        <outputdir> is the directory to save png files to 
        """
        print """\nExample:
        produce month, year and all dbploteventrate plots from the AVO Catalog for each volcano in the volcanoes_avo database \n
        %s /Seis/Kiska4/picks/Total/Total volcanoes_avo /usr/local/mosaic/AVO/avoseis/counts
        """ % (sys.argv[0])

def main(argv=None):
        try:
                opts, args = getopt.getopt(argv, 'vh')
                if len(args)<3:
                        usage()
                        sys.exit(2)
                catalogpath = args[0]
                dbplacespath = args[1]
                outdir = args[2]
        except getopt.GetoptError,e:
                print e
                usage()
                sys.exit(2)

        verbose = False
        savepercentiles = True

        for o, a in opts:
                if o == "-v":
                        verbose = True
                elif o in ("-h", "--help"):
                        usage()
                        sys.exit()
                else:
                        assert False, "unhandled option"

        if verbose:
                print "catalogpath = " + catalogpath
                print "dbplacespath = " + dbplacespath
                print "outdir = " + outdir

        # time now
        datetimenow = datetime.datetime.now() # datetime
        #epochnow = time.mktime(datetimenow.timetuple()) # epoch
        epochnow = datascope.stock.now()
        #datenumnow = mpl.dates.epoch2num(epochnow) # datenumber
        secsperday = 60 * 60 * 24
	daysperyear = 365
	dayspermonth = 30
        epoch1989 = 599616000

	# epoch at start and end today
	epochtodaystart = (epochnow // secsperday) * secsperday
	epochtodayend = epochtodaystart + secsperday
	
	# datenum start and end time
	enum = mpl.dates.epoch2num(epochtodayend)
	epoch1yearago = epochtodaystart - (secsperday * daysperyear)
	epoch1monthago = epochtodaystart - (secsperday * dayspermonth)
	snum1yearago = mpl.dates.epoch2num(epoch1yearago)
	snum1monthago = mpl.dates.epoch2num(epoch1monthago)
	
        dictplaces = modgiseis.readplacesdb(dbplacespath)
        place = dictplaces['place']
        lat = dictplaces['lat']
        lon = dictplaces['lon']
        radius = dictplaces['radius']
	n = place.__len__()
	if n > 0:
		print "- number of places = {}".format(n)
	
		for c in range(n):
	
			# LAST MONTH
			monthfile = outdir + "/" + place[c] + "_month.png"
			if os.path.exists(monthfile):
				os.remove(monthfile)
			subset_expr = "time > %f && deg2km(distance(lat, lon, %s, %s))<%s" % (epoch1monthago, lat[c], lon[c], radius[c])
			dbploteventrate.main([catalogpath, monthfile, subset_expr, snum1monthago, enum])
	
	
			# LAST YEAR
			yearfile = outdir + "/" + place[c] + "_year.png"
			if os.path.exists(yearfile):
				os.remove(yearfile)
			subset_expr = "time > %f && deg2km(distance(lat, lon, %s, %s))<%s" % (epoch1yearago, lat[c], lon[c], radius[c])
			dbploteventrate.main([catalogpath, yearfile, subset_expr, snum1yearago, enum])
	
			# TOTAL
			totalfile = outdir + "/" + place[c] + "_total.png"
			if os.path.exists(totalfile):
				os.remove(totalfile)
			subset_expr = "time > %f && deg2km(distance(lat, lon, %s, %s))<%s" % (epoch1989, lat[c], lon[c], radius[c])
			dbploteventrate.main([catalogpath, totalfile, subset_expr])


        ############

        print "Done.\n"


if __name__ == "__main__":
        if len(sys.argv) > 1:
                main(sys.argv[1:])
        else:
                usage()


