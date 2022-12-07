#!/usr/bin/env python
import matplotlib.pyplot as plt
import pandas as pd
from random import gauss
from random import seed
from pandas import Series
import os, sys
if len(sys.argv) < 2:
    print("Usage: %s DEPTH RADIUS\n" % sys.argv[0])
z = sys.argv[1]
a = sys.argv[2]
filename = "dipole_%s_%s.csv" % (z, a)

colnames = ['x', 'M']
#data = pd.read_csv(filename, names=colnames, header=None)
data = pd.read_csv(filename)
data.columns = data.columns.str.replace(' ', '')

# create white noise series
seed(1)
magamp = 0.05 # nT at 100 Hz sampling rate
series = [gauss(0.0, magamp) for i in range(len(data.index))]
series = Series(series)

# should add environmental noise too

plt.scatter(data['x'], data['M'] + series)
plt.xlabel('Distance along profile (m)')
plt.ylabel('Magnetic anomaly (nT)')
plt.title('Buried sphere, depth %s, radius %s' % (z, a))
pngfile = filename.replace("csv", "png")
plt.savefig(pngfile)
os.system("open %s" % pngfile)
