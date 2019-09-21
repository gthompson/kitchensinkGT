#!/raid/apps/OBSPY/bin/python
from obspy.core import read
st = read('http://examples.obspy.org/RJOB_061005_072159.ehz.new')
print(st)
len(st)
tr = st[0]  
print(tr)
print(tr.stats)
tr.stats.station
tr.data
len(tr)
st.plot(outfile='example_plot.png')
