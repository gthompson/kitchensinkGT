from obspy.core import UTCDateTime, read, stream
from obspy.signal import seisSim, invsim, detrend
import matplotlib.pyplot as plt
import os
import numpy as np
from streamplot import streamplot
RESPDIR = "/home/t/thompsong/newtonhome/seismo/CAL"

def instrumentcorrect(rawstream, boolPlot=0):
    # the corrected stream to be returned
    correctedstream = stream.Stream()
    
    #fl1=0.00588
    #fl2=0.00625
    fl1=0.3
    fl2=0.5
    fl3=20.
    fl4=33.
    
    # loop over all stream objects
    n = len(rawstream)
    for i in range(n):
        print("Channel %d:" % i)
        channelstream = stream.Stream(traces=[rawstream[i].copy()])
        print(channelstream[0])
        # This works for 2000 data except for pressure sensor and MBMH (which isn't in dataless seed from Silvio)
            if channelstream[0].stats.channel == 'SBZ':
                    channelstream[0].stats.channel = 'BHZ'
            elif channelstream[0].stats.channel == 'SBE':
                    channelstream[0].stats.channel = 'BHE'
            elif channelstream[0].stats.channel == 'SBN':
                    channelstream[0].stats.channel = 'BHN'
            elif channelstream[0].stats.channel == 'S Z':
                    channelstream[0].stats.channel = 'SHZ'
            elif channelstream[0].stats.channel == 'S E':
                    channelstream[0].stats.channel = 'SHE'
            elif channelstream[0].stats.channel == 'S N':
                    channelstream[0].stats.channel = 'SHN'
            elif channelstream[0].stats.channel == 'BH':
                    channelstream[0].stats.channel = 'BH' + channelstream[0].stats.location
            elif channelstream[0].stats.channel == 'SH':
                    channelstream[0].stats.channel = 'SH' + channelstream[0].stats.location
        print(channelstream[0].stats)
    
        respf = RESPDIR + '/' + 'RESP.MN.' + channelstream[0].stats.station + '..' + channelstream[0].stats.channel
        if os.path.exists(respf):
            print(respf + " found")
        else:
            print(respf + " not found")
            continue

        # TRACE 1 (channelstream[0]) is raw seismogram
    
        # TRACE 2 (channelstream[1]) is detrended seismogram
        # remove last second (there is a spike in last second of most files it seems), then detrend
        channelstream.append(channelstream[0].copy())
        fsamp = channelstream[1].stats.sampling_rate
        y = channelstream[1].data
        y[len(y)-int(fsamp/2):len(y)] = np.median(y)
        channelstream[1].data = y
        channelstream[1].detrend()
    
        # TRACE 3 (channelstream[2]) is corrected seismogram    
        channelstream.append(channelstream[1].copy())
        print("Correcting for instrument response")
        seedresp = {'filename':respf,'date':channelstream[2].stats.starttime, 'units':'VEL'}
        #channelstream[2].data = seisSim(channelstream[2].data, fsamp, paz_remove=None, remove_sensitivity=False, pre_filt=(fl1,fl2,fl3,fl4), seedresp=seedresp)
        channelstream[2].data = seisSim(channelstream[2].data, fsamp, paz_remove=None, pre_filt=(fl1,fl2,fl3,fl4), seedresp=seedresp)

        # PLOT THE TRACE VECTOR
        if boolPlot:
            streamplot(channelstream)

        # ADD THE CORRECTED TRACE TO THE CORRECTEDSTREAM OBJECT
        correctedstream.append(channelstream[2])

    return correctedstream

if __name__ == "__main__":
    #rawstream = read('/raid/data/seisan/WAV/MVOE_/2000/01/2000-01-14-1237-35S.MVO___019')
    rawstream = read('/raid/data/seisan/WAV/MVOE_/2001/07/2001-07-01-1216-04S.MVO___018')
    correctedstream = instrumentcorrect(rawstream, 1)
    print("Got %d corrected traces from %d raw traces" % (len(correctedstream), len(rawstream)))
    streamplot(correctedstream)

