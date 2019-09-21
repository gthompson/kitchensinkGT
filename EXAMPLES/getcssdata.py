# getting waveform data
from obspy.core import read
dbpath='/opt/antelope/data/db/demo/demo.wfdisc'
st = read(dbpath, format='CSS')
st.plot()

# getting catalog data
c = CSS_Catalog.CSS_Catalog()
c.read('/opt/antelope/data/db/demo/demo')
print c.events
print c.events[0].origins[0]

