from obspy.core import *
import sys
for file in sys.argv[1:]:
    st = read(file)
    for tr in st:
	outfile = "%s.%s.%s" % (tr.stats.station, tr.stats.channel, tr.stats.starttime)
	tr.write("outfile + ".sac", format="SAC")
