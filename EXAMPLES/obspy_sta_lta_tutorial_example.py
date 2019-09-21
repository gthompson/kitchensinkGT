# -*- coding: utf-8 -*-
"""
Geophysical Time Series Analysis Week 14
Examples of how to do single-station & multi-station STA/LTA
triggering in Python with ObsPy.
Includes examples on Pavlof 2007 explosion-quakes, and Redoubt
2009 swarm.

Glenn Thompson April 2016
"""

# 1. Load trace (waveform) which has sampling_rate = 200 Hz
from obspy.core import read
st = read("https://examples.obspy.org/ev0_6.a01.gse2")
st = st.select(component="Z")
tr = st[0]

# 2. Show information about this trace
tr.stats
df = tr.stats.sampling_rate

# 3. Plot the time series!
tr.plot(type="relative")

# 4. Compute classic_sta_lta with STA=1000 samples (5s) and LTA=2000 samples (10s)
from obspy.signal.trigger import classic_sta_lta
help(classic_sta_lta)
staltaratio = classic_sta_lta(tr.data, int(5*df), int(10*df))

# 5. Plot the trace with corresponding STA/LTA ratio
from obspy.signal.trigger import plot_trigger
help(plot_trigger)
plot_trigger(tr, staltaratio, 1.5, 0.5)

# 6. Repeat 4 & 5, but with Z-detect
from obspy.signal.trigger import z_detect
help(z_detect)
staltaratio = z_detect(tr.data, int(df*10))
plot_trigger(tr, staltaratio, -0.4, -0.3)

# 7. Repeat 4 & 5, but with recursive_sta_lta
from obspy.signal.trigger import recursive_sta_lta
staltaratio = recursive_sta_lta(tr.data, int(5 * df), int(10 * df))
help('recursive_sta_lta')
plot_trigger(tr, staltaratio, 1.2, 0.5)

# 8. Repeat 4 & 5, but with carl_sta_trig (Earthworm)
from obspy.signal.trigger import carl_sta_trig
help('carl_sta_trig')
staltaratio = carl_sta_trig(tr.data, int(df * 5), int(10 * df), 0.8, 0.8)
plot_trigger(tr, staltaratio, 20.0, -20.0)

# 9. Repeat 4 & 5, but with delayed_sta_lta
from obspy.signal.trigger import delayed_sta_lta
staltaratio = delayed_sta_lta(tr.data, int(5*df), int(10*df))
plot_trigger(tr, staltaratio, 5, 10)

#####################################
# 10. Let's try to optimize the STA & LTA settings for this data
# To do this I wrote a function "tune_sta_lta". Let's view the help for this
# function. First we have to find and import the module.
import sys
sys.path.append('/Users/glennthompson/Dropbox/scratch_matlab')
import tune_sta_lta as tsl
help(tsl.tune_sta_lta)

# read a seismogram into a trace object
from obspy.core import read
st = read("https://examples.obspy.org/ev0_6.a01.gse2")
st = st.select(component="Z")
tr = st[0]

# plot time series
tr.plot()

# plot spectrogram
tr.spectrogram()

# run STA/LTA 100 times to find best settings
algorithm = 'classic_sta_lta'
TSIGNAL_START = 30.0
TSIGNAL_END = 40.0
NTRIES=100
tsl.tune_sta_lta(tr, algorithm, TSIGNAL_START, TSIGNAL_END, NTRIES)

