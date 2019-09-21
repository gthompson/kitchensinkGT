#!/raid/apps/OBSPY/bin/python
# Simple stream object plotter
# Glenn Thompson 2013/01/07

import numpy as np
import matplotlib.pyplot as plt

def streamplot(st, bottomlabel='', ylabels=[]):

	# close all plots
	plt.figure()

	# start time as a Unix epoch
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
			# for the bottom subplot, also add an xlabel with start time
			if bottomlabel=='':
				plt.xlabel("Starting at %s" % (st[0].stats.starttime) )
			else:
				plt.xlabel(bottomlabel)

		# default ylabel is station.channel
		if ylabels==[]:
			plt.ylabel(st[i].stats.station + "." + st[i].stats.channel, rotation=0)
		else:
			plt.ylabel(ylabels[i])

		# explicitly give the maximum amplitude and offset(median)
		plt.text(0, 1, "max=%.1e offset=%.1e" % (np.max(np.abs(y)), offset),
        		horizontalalignment='left',
        		verticalalignment='top',transform=axh[i].transAxes)

	# change all font sizes
	plt.rcParams.update({'font.size': 8})

	# show the figure
	plt.show()

