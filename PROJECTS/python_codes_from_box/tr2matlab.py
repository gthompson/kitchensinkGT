from scipy.io import savemat

import obspy


st = obspy.read("https://examples.obspy.org/BW.BGLD..EH.D.2010.037")
for i, tr in enumerate(st):
    mdict = {k: str(v) for k, v in tr.stats.iteritems()}
    mdict['data'] = tr.data
    #mdict['sampling_rate'] = tr.stats.sampling_rate
    #mdict['network'] = tr.stats.network
    #mdict['station'] = tr.stats.station
    #mdict['location'] = tr.stats.location
    #mdict['channel'] = tr.stats.channel
    #mdict['starttime'] = tr.stats.starttime
    #mdict['endtime'] = tr.stats.endtime
    #mdict['delta'] = tr.stats.delta
    #mdict['npts'] = tr.stats.npts
    #mdict['calib'] = tr.stats.calib
    #mdict['format'] = tr.stats._format

    savemat("data-%d.mat" % i, mdict)
