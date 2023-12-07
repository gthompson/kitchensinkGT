from obspy.core import read
from datetime.datetime import UTCDateTime
st = read('/raid/data/rockets/wfdata/rocketevents/20160901_SpaceXplosion/mseed/FL*')
starttime = UTCDateTime("2016-09-01T13:07:00")
endtime = UTCDateTime("2016-09-01T13:10:00")
st.trim(starttime, endtime)
st.write("explosion.mseed", format="MSEED")

