from obspy.core import UTCDateTime
from obspy.earthworm import Client
client = Client("pubavo1.wr.usgs.gov", 16022)
response = client.availability("AV", "AKV", channel="EHZ")
print(response)  
t = response[0][4]
tstart = UTCDateTime(2013,9,20,0o5,30,0)
tend = UTCDateTime(2013,9,20,0o5,40,0)
if tstart >= response[0][4] and tend <= response[0][5]:
    st = client.getWaveform('AV', 'AKV', '', 'EHZ', tstart, tend)
    st.plot()
 
