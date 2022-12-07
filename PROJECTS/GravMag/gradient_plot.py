#!/usr/bin/env python
import matplotlib.pyplot as plt
import pandas as pd
from random import gauss
from random import seed
from pandas import Series
import os, sys
if len(sys.argv) < 3:
    print("Usage: %s DEPTH1 DEPTH2 RADIUS\n" % sys.argv[0])
z1 = sys.argv[1]
z2 = sys.argv[2]
a = sys.argv[3]
filename1 = "dipole_%s_%s.csv" % (z1, a)
filename2 = "dipole_%s_%s.csv" % (z2, a)

#colnames = ['x', 'M']
#data = pd.read_csv(filename, names=colnames, header=None)
data1 = pd.read_csv(filename1)
data1.columns = data1.columns.str.replace(' ', '')
data2 = pd.read_csv(filename2)
data2.columns = data2.columns.str.replace(' ', '')

# create white noise series
seed(1)
magamp = 0.05 # nT at 100 Hz sampling rate
series1 = [gauss(0.0, magamp) for i in range(len(data1.index))]
series1 = Series(series1)
series2 = [gauss(0.0, magamp) for i in range(len(data2.index))]
series2 = Series(series2)

M1 = data1['M'] + series1
M2 = data2['M'] + series2
MZdiff = M1 - M2

f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(data1['x'], data1['M']-data2['M'], label='Vertical')
ax1.plot(data1['x'], data1['M'].diff(periods=10), label='Horizontal')
#ax1.set_xlabel('Distance along profile (m)')
ax1.set_ylabel('Magnetic gradient\n anomaly (nT)')
ax1.legend()
ax2.scatter(data1['x'], MZdiff, 5, label='Vertical')
ax2.scatter(data1['x'], M1.diff(periods=10), 5, label='Horizontal')
ax2.set_xlabel('Distance along profile (m)')
ax2.set_ylabel('Magnetic gradient\n anomaly (nT)')
ax2.legend()
ax1.set_title('Buried sphere, h1 = %s m, h2 = %s m, a = %s m\nTheoretical' % (z1, z2, a))
ax2.set_title('Sampled every 10-cm with 0.05 nT instrument noise')
ax1.set_xlim(-50, 50)
ax2.set_xlim(-50, 50)
pngfile = "gradient_%s_%s_%s.png" % (z1, z2, a)
f.savefig(pngfile)
os.system("open %s" % pngfile)



