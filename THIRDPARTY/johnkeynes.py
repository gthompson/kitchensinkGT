from obspy import read
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime

kwargs = dict(water_level=60, pre_filt=(1, 3, 20, 30))
resp="/home/hd/Dropbox/projects/current/RESP.FABD.FA.6N.HNE"
wavfile="2007221075749.60.FABD.HNE_6N.gz"

plt.figure(1)
plt.subplot(311)
plt.title('Uncorrected')
tr = read(wavfile)
plt.plot(tr[0].data)

plt.subplot(312)
plt.title('Corrected')
tr = read(data)
tr.simulate(seedresp={"filename": resp,'date': UTCDateTime("2007-10-22T07:57:49.000"),"units": "ACC"},**kwargs)
#tr.filter(type="highpass", freq=0.1)
plt.plot(tr[0].data)
plt.subplot(313)
plt.title('Corrected via RESP file for [DIS,VEL+integrate,ACC+2xintegrate]+highpass filter')
tr = read(data)
tr.simulate(seedresp={"filename": resp,'date': UTCDateTime("2007-10-22T07:57:49.000"),"units": "DIS"},**kwargs)
#tr.filter(type="highpass", freq=0.1)
plt.plot(tr[0].data)
tr = read(data)
tr.simulate(seedresp={"filename": resp,'date': UTCDateTime("2007-10-22T07:57:49.000"),"units": "VEL"},**kwargs)
tr.integrate()
#tr.filter(type="highpass", freq=0.1)
plt.plot(tr[0].data)
tr = read(data)
tr.simulate(seedresp={"filename": resp,'date': UTCDateTime("2007-10-22T07:57:49.000"), "units": "ACC"},**kwargs)
tr.integrate()
tr.integrate()
tr.filter(type="highpass", freq=0.1)
plt.plot(tr[0].data)
plt.savefig("plot.pdf",format='pdf')
plt.savefig("plot.eps",format='eps')
plt.savefig("plot.png",format='png')    
