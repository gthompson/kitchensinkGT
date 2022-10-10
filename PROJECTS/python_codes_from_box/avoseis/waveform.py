#! /usr/bin/env python
""" waveform class, Glenn Thompson 2010/05/01 """

def hanningcosine(data, fraction=1/3):
    """ create a cosine tapered hanning window with a flat central section """
    from numpy import round 
    numsamplesdata = len(data)
    numsamplestaper = fraction * round(numsamplesdata)
    from scipy import signal
    h = signal.hanning(numsamplestaper * 2)
    h1 = h[:numsamplestaper]
    h2 = h[-numsamplestaper]
    win = [h1, [1.0]*numsamplesdata, h2]
    return win

class waveform(object):
    
    def __init__(self, estart = 0, samprate = 100.0, sta = 'DUMM', chan = 'EHZ', net = 'AV', data=[], seismogramtype = 'raw'):
        """ initialise an empty waveform object """
	self.estart = estart
	self.samprate = samprate
	self.sta = sta
	self.chan = chan
	self.net = net
	self.data = data
	self.seismogramtype = seismogramtype

    def __str__(self):
        """ print out the attributes of a waveform object """
	print 'sta: %s' % self.sta
	print 'chan: %s' % self.chan
	print 'net: %f' % self.net
        print 'start: %f' % self.estart
	print 'end: %f' % self.estart + len(self.data)/self.samprate
	print 'samprate: %f' % self.samprate
	print 'number of samples: %d', len(self.data)
            
    def export_sam(self):
	"""
	take waveform data and compute sam
	waveform data can have 'seismogramtype' of:
		'raw' -> rsam
		'velocity' seismogram -> vsam
		'displacement' seismogram -> dsam
		energy -> esam (use formula from Brando)
	"""
	name = {'raw' -> 'rsam', 'velocity' -> 'vsam', 'displacement' -> 'dsam', 'energy' -> 'esam'}
	thisname = name{self.seismogramtype}
	print 'Computing %s' % thisname
	numminutes = (len(self.data)/self.samprate)/60
	c = 0
	for i in range(numminutes):
		wdatasubset = self.data[c:c+self.samprate*60-1]
		if thisname != 'esam':
			metric[i] = nanmedian(wdatasubset)
		else:
			metric[i] = nansum(wdatasubset)
		c += self.samprate*60
	import sam
	eend = self.estart + len(self.data)/self.samprate
	samobject = sam.sam(thisname, self.estart, eend, self.sta, self.chan, '', metric)
	
		
    def reconstitute(self, deconvolve_flag=FALSE, integrate_flag = FALSE):
	"""
	reconstitute waveform
	first filter it robustly between 0.8 and 15 Hz
	then apply calibration correction (or deconvolve if deconvolve_flag)
	"""
	self.detrend()
	self.filter(0.8, 15, 4)
	if deconvolve_flag:
		self.deconvolve()
	else:
		self.calibrate()
	if integrate_flag:
		self.integrate()
		self.seismogramtype = 'displacement'
	else:
		self.seismogramtype = 'velocity'

   def filter(self, lowcut=0.8, highcut=15.0, poles=4):
	"""
	filter
	# create d = [wdata-reversed][wdata][wdata-reversed]
	# apply 50% hanning window
	# apply strong bandpass filter 0.8 - 15 Hz, 2 way
	# wdata = extract middle 3rd	
	"""
	wdata = self.data
	numsamples = len(wdata)
	wdata3 = [fliplr(wdata), wdata, fliplr(wdata)]
	win = hanningcosine(wdata3)
	wdata3 = wdata3 .* win
	# create bandpass filter lowcut to highcut
	# apply both ways on wdata3
	# extract middle 3rd
	wdata = wdata3[numsamples + 1: numsamples * 2]
	return wdata
	



	
		
		