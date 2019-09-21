from obspy.core import UTCDateTime
from obspy.earthworm import Client
client = Client("pele.ess.washington.edu", 16017)
response = client.availability("UW", "TUCA", channel="BHZ")
print(response)  
t = response[0][4]
st = client.getWaveform('UW', 'TUCA', '', 'BH*', t + 100, t + 130)
st.plot() 
