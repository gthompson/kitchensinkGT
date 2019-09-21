# Create a total database from a list of daily databases
# Glenn Thompson, 2012-08-10

import glob, shutil
sys.path.append(os.environ['ANTELOPE'] + "/data/python")
import antelope.datascope as datascope

if len(sys.argv) < 3:
	# stop program and print an error message
	sys.exit("Usage: %s /path/pattern/to/source/databases /path/to/target/total/database [/path/to/current/db]\ne.g.\n\t%s events/optimised/events_????_??_?? events/optimised/Total events/optimised/events" % (sys.argv[0], sys.argv[0]))

sourcedbpattern = sys.argv[1]
targetdb = sys.argv[2]
datascope.dbcreate(targetdb, "css3.0")
for table in ['origin', 'event', 'assoc', 'arrival', 'netmag', 'stamag']:
	destination = open("%s.%s" % (targetdb, table), 'wb')
	print "SOURCE = " + sourcedbpattern + "." + table
	print "TARGET = " + targetdb + "." + table + '\n'
	for filename in glob.glob("%s.%s" % (sourcedbpattern, table)):
		print filename
		if os.path.islink(filename):
			print "%s is a link - skipping" % filename
		else:
			print "will concatenate %s " % filename
			shutil.copyfileobj(open(filename, 'rb'), destination)
	if len(sys.argv) == 4:
		if os.path.isfile(sys.argv[3]):
			shutil.copyfileobj(open("%s.%s" % (sys.argv[3],table), 'rb'), destination)
	destination.close()
		
	
