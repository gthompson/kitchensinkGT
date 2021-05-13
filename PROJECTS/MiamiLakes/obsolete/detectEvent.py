from obspy.core import read
import libseisGT
import matplotlib.pyplot as plt
from obspy.signal.trigger import classic_sta_lta, z_detect, plot_trigger, trigger_onset
import glob
mseedfiles = glob.glob("eventMiniseed/Event*.CORRECTED.mseed")
pretrig = 30
posttrig = 10
for mseed in mseedfiles:

    st = read(mseed)
    st = libseisGT.detectEvent()
    isofmt = st[0].stats.starttime.isoformat()
    seisanmseed = isofmt[0:10] + "-" + isofmt[11:13] + isofmt[14:16] + "-" + isofmt[17:19] + 'S.MILAK__%03d' % (len(st)+1)
    print('- writing %s' % seisanmseed)
    pngfile = os.path.join(TOPDIR, 'events', 'plots', 'detected', seisanmseed + '.png')
    st.plot(outfile = pngfile, equal_scale = False)
    st.write(os.path.join(TOPDIR, 'events', 'mseed', 'detected', seisanmseed), format='MSEED' )
print('Now run dirf and autoreg to put these into the MILAK Seisan database')
