#!/raid/apps/OBSPY/bin/python
# Show de-medianed seismograms in the manner of a record section like mulplt
# Glenn Thompson 2013/01/06

import numpy as np
import matplotlib.pyplot as plt
from obspy import read
import datetime, sys


def main():

	# read in the wav file name from the command line
	wavfile = sys.argv[1]

	# close all plots
	plt.close('all')

	# read the wav file as a stream object
	st = read(wavfile)
	print st

	#fh = plt.figure()

	# Only do this if want to use a common x-axis
	#ax = plt.subplot(n, 1, 1)

	# start time of the record section
	startepoch = st[0].stats.starttime.timestamp

	# create empty set of subplot handles - without this any change to one affects all
	axh = []

	# loop over all stream objects
	n = len(st)
	for i in range(n):
		# add new axes handle for new subplot
		#axh.append(plt.subplot(n, 1, i+1, sharex=ax))
		axh.append(plt.subplot(n, 1, i+1))

		# time vector, t, in seconds since start of record section
    		t = np.linspace(st[i].stats.starttime.timestamp - startepoch,
                    	st[i].stats.endtime.timestamp - startepoch,
                    	st[i].stats.npts)

		# We could detrend, but in case of spikes, subtracting the median may be better
		#st[i].detrend()
		offset = np.median(st[i].data)
		y = st[i].data - offset

		# PLOT THE DATA
		axh[i].plot(t, y)

		# remove yticks because we will add text showing max and offset values
		axh[i].yaxis.set_ticks([])

		# remove xticklabels for all but the bottom subplot
		if i < n-1:
			axh[i].xaxis.set_ticklabels([])
		else:
			# for the bottom subplot, also add an xlabel with wavfilename and start time
			plt.xlabel("WAV file %s\n Starting at %s" % (wavfile, st[0].stats.starttime) )

		# ylabel is station.channel
		plt.ylabel(st[i].stats.station + "." + st[i].stats.channel, rotation=0)

		# explicitly give the maximum amplitude and offset(median)
		plt.text(0, 1, "max=%.1e offset=%.1e" % (np.max(np.abs(y)), offset),
        		horizontalalignment='left',
        		verticalalignment='top',transform=axh[i].transAxes)

		#st[i].plot(fig=fh, color='black', tick_format='%I:%M:%S.%s', starttime=st[i].stats.starttime, endtime=st[i].stats.starttime+20)

	# change all font sizes
	plt.rcParams.update({'font.size': 8})

	# show the figure
	plt.show()

	#st[0].spectrogram(log=True, title=st[0].stats.station + " " + str(st[0].stats.starttime))


if __name__ == "__main__":
    main()
