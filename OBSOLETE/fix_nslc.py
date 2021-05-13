#!/raid/apps/anaconda2/bin/python
''' an old code that changed BHZ to HHZ, BDF to BHF etc.
and changed location code from 00 to 10 at 2016/5/1 and 20 at 2017/1/1.
presumably the location changes were for the Beach House station '''
from obspy import read
import os.path
import numpy as np
import datetime
import sys
network = 'FL'

file=open('wffiles.txt','r')
lines = file.readlines()
#os.system('cat ' + filenrlis)
for thisfile in lines:
        print thisfile

        st = read(thisfile[:-1])
        #try:
        #        st = read(thisfile)
        #except:
        #        print "Cannot read " + thisfile
        #        continue

        for tr in st:


          # fix the network, channel and location
                print " "
                print tr.stats
                tr.stats['network']=network
                sta = tr.stats['station'].strip()
                chan = tr.stats['channel'].strip()

		# fix B channels
                if chan[0:1]=='B':
                        chan='H'+chan[1:3]

		# fix time periods
		dt1 = datetime.datetime(2016,5,1,0,0,0);		
		dt2 = datetime.datetime(2017,1,1,0,0,0);
                starttime = tr.stats['starttime']
		location = '00'
		if starttime > dt1:
			location = '10'
		if starttime > dt2:
			location = '20'

                jjj = str(starttime.julday).zfill(3)
                yyyy = str(starttime.year).zfill(4)
                mm = str(starttime.month).zfill(2)
                dd = str(starttime.day).zfill(2)
                hh = str(starttime.hour).zfill(2)
                ii = str(starttime.minute).zfill(2)
                ss = str(starttime.second).zfill(2)

		# now fix tr.stats
                tr.stats['network'] = network
                tr.stats['channel'] = chan
                tr.stats['location'] = location

		print "Changed to:"
		print tr.stats

                # Write Miniseed files
                antelopedir = yyyy + "/" + jjj
                if not os.path.exists(antelopedir):
                        os.makedirs(antelopedir)
                antelopefile = sta + "." + chan + "." + yyyy + mm + dd + "T" + hh + ii + ss
                antelopefullpath = antelopedir + "/" + antelopefile
                print antelopefullpath
                tr.write(antelopefullpath, format="MSEED")

