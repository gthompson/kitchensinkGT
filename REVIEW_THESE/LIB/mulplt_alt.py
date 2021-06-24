#!/raid/apps/OBSPY/bin/python
# Show detrended seismograms in the manner of a record section
import numpy as np
import matplotlib.pyplot as plt
from obspy import read
import datetime

# read in the wav file name from the command line
import sys
wavfile = sys.argv[1]

# close all plots
plt.close('all')

# read the wav file as a stream object
st = read(wavfile)
print(st)

#fh = plt.figure()

# Only do this if want to use a common x-axis
#ax = plt.subplot(n, 1, 1)


startepoch = st[0].stats.starttime.timestamp

# create empty set of subplot handles - without this any change to one affects all
axh = []

# loop over all stream objects
n = len(st)
for i in range(n):
    #axh.append(plt.subplot(n, 1, i+1, sharex=ax))
    axh.append(plt.subplot(n, 1, i+1))
    #st[i].detrend()
        t = np.linspace(st[i].stats.starttime.timestamp - startepoch,
                        st[i].stats.endtime.timestamp - startepoch,
                        st[i].stats.npts)
    offset = np.median(st[i].data)
    axh[i].plot(t, st[i].data - offset)

    axh[i].yaxis.set_ticks([])
    if i < n-1:
        axh[i].xaxis.set_ticklabels([])
    else:
        plt.xlabel("WAV file %s\n Starting at %s" % (wavfile, st[0].stats.starttime) )
    plt.ylabel(st[i].stats.station + "." + st[i].stats.channel, rotation=0)
    plt.text(0, 1, "max=%.1e offset=%.1e" % (np.max(st[i].data), offset),
            horizontalalignment='left',
            verticalalignment='top',transform=axh[i].transAxes)

    #st[i].plot(fig=fh, color='black', tick_format='%I:%M:%S.%s', starttime=st[i].stats.starttime, endtime=st[i].stats.starttime+20)


plt.rcParams.update({'font.size': 8})


plt.show()
#st[0].spectrogram(log=True, title=st[0].stats.station + " " + str(st[0].stats.starttime))

