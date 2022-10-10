#! /usr/bin/env python
""" sam class, Glenn Thompson 2010/05/01 """
""" Goal is to be able to have an M*N matrix of sam objects - and then plot them accordingly 
So need to init an M*N matrix of sam objects 
Then load an M*N matrix
Filter them
Plot them
"""

class sam(object):
	
    import datetime
    now = datetime.datetime.utcnow()
    
    def __init__(self, name='rsam', estart = 0, eend = 360, sta = 'DUMM', chan = 'EHZ', filepattern = 'F:/RSAM/STACHANYYYY.DAT', data=[]):
        """ initialise an empty sam object """
        self.name = name
	self.estart = estart
	self.eend = eend
	self.sta = sta
	self.chan = chan
	self.filepattern = filepattern
	self.data = data

    def __str__(self):
        """ print out the attributes of a sam object """
        print '%s' % self.name,' %12d' % self.estart,' %12d' % self.eend,' %s' % self.sta, ' %s' % self.chan, ' %s' % self.filepattern
            
    def import_samfile(self):
	"""
	open file
	work out number of samples to get
	position pointer
	read samples
	close file
	return data
	Examine:
	numpy.fromfile (file=, dtype=int, count=-1, sep=”)
Return a 1-d array of data type, dtype, from a file (open file object or string with
the name of a file to read). The file will be read in binary mode if sep is the
empty string. Otherwise, the file will be read in text mode with sep providing
the separator string between the entries. If count is -1, then the size will be
determined from the file, otherwise, up to count items will be read from the
file. If fewer than count items are read, then a RunTimeWarning is issued
indicating the number of items read.
	
	"""
	
	
    def export_samfile(self):
	"""
	open file
	work out number of samples to write
	position pointer
	write samples
	close file
	"""	
	
    def import_samfile_wrapper(self):
	"""
	loop over year files:
		import_samfile
		extend data list
	return data
	"""	
	
    def create_samfile(self):
	"""
	open file for output
	write 366 days worth of 0
	end
	"""
	
    def downsample_sam(self, factor):
	"""
	take every factor sample
	"""	

	
    def plot_sam(self):
	"""
	plot_date sam
	"""	

 