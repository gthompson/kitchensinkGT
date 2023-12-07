import matplotlib.pyplot as plt  
from obspy.core import read
import numpy as np
import mlpy

stream = read("http://examples.obspy.org/a02i.2008.240.mseed")
trace = stream[0]

omega0 = 8
spec, scale = mlpy.cwt(trace.data, dt=trace.stats.delta, dj=0.05,
                        wf='morlet', p=omega0, extmethod='none',
                        extlength='powerof2')
freq = (omega0 + np.sqrt(2.0 + omega0**2)) / (4*np.pi * scale[1:])

t = np.arange(trace.stats.npts) / trace.stats.sampling_rate
plt.imshow(np.abs(spec), extent=(t[0], t[-1], freq[-1], freq[0]))
plt.xlabel('Time [s]')
plt.ylabel('Frequency [Hz]')
plt.show()
