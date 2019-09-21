from obspy.core import read
st = read("/Volumes/GoogleDrive/My Drive/data/KSC/explosion.mseed")
st[0].stats
#st.plot()
#stime = UTCDateTime("2002-10-23T11:27:00")
#etime = UTCDateTime("2002-10-23T11:29:00")
#st.plot(starttime=stime, endtime=etime)
