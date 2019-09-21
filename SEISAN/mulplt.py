#!/raid/apps/OBSPY/bin/python
# Show de-medianed seismograms in the manner of a record section like mulplt
#from obspy import read
import sys
import StreamGT
import Seisan_Catalog
#from streamplot import streamplot

def main():

	# read in the wav file name from the command line
	if len(sys.argv) > 1:
		path = sys.argv[1]
	else:
		print "*** TEST MODE ***\nNo wavfile given on command line, so using a default"
		path = '/raid/data/seisan/WAV/MVOE_/2000/01/2000-01-14-1237-35S.MVO___019'

	# read the wav file as a stream object
        w = Wavfile(path)
        w.read()
        w.plot()
        w.mulplt()

	#st = read(path)
	#print st
	#streamplot(st)

if __name__ == "__main__":
	main()
