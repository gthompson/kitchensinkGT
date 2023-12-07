from obspy.iris import Client
from obspy import UTCDateTime
client = Client()

# get waveform data
t = UTCDateTime("2010-02-27T06:45:00.000")
st = client.getWaveform("IU", "ANMO", "00", "BHZ", t, t + 60 * 60)
st.plot()  

# save waveform data
t1 = UTCDateTime("2010-02-27T06:30:00.000")
t2 = UTCDateTime("2010-02-27T07:30:00.000")
client.saveWaveform('IU.ANMO.00.BHZ.mseed', 'IU', 'ANMO',
                    '00', 'BHZ', t1, t2) 

# save response data
#t = UTCDateTime(2009, 1, 1)
#client.saveResponse('resp.txt', 'IU', 'ANMO', '', '*',
#                    t, t + 1, format="RESP") 

# get events
starttime = UTCDateTime("2011-04-01")
endtime = UTCDateTime("2011-04-15")
cat = client.getEvents(starttime=starttime, endtime=endtime,
                       minmag=6.7)
print(cat)  
cat.plot()  
