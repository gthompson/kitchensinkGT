#! /usr/bin/env python
""" sam class, Glenn Thompson 2010/05/01
    Goal is to be able to have an M*N matrix of sam objects
    and then plot them accordingly 
    So need to init an M*N matrix of sam objects 
    Then load an M*N matrix
    Filter them
    Plot them
"""
from obspy.core import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt

class sam(object):
    
    def __init__(self, metric, starttime, endtime, 
        station, channel, data):
        """ initialise an empty sam object """
        self.metric = metric
        self.starttime = starttime
        self.endtime   = endtime
        self.station = station
        self.channel = channel
        self.data = data

    def __str__(self):
        """ print out the attributes of a sam object """
        str =  'sam object:\n'
        str += '\tmetric:    %s\n'   % self.metric
        str += '\tstarttime: %r\n'   % self.starttime
        str += '\tendtime:   %r\n'   % self.endtime
        str += '\tstation:   %s\n'   % self.station
        str += '\tchannel:   %s\n'   % self.channel
        str += '\tdata:      %r\n'   % self.data
        return str

    def import_samfile(self):
        """
        open file
        work out number of samples to get
        position pointer
        read samples
        close file
        return data
        Examine:
        numpy.fromfile(file=, dtype=int, count=-1, sep='')
        Return a 1-d array of data type, dtype, 
        from a file (open file object or string with
        the name of a file to read). The file will be 
        read in binary mode if sep is the
        empty string. Otherwise, the file will be read 
        in text mode with sep providing
        the separator string between the entries. If 
        count is -1, then the size will be
        determined from the file, otherwise, up to count 
        items will be read from the
        file. If fewer than count items are read, then 
        a RunTimeWarning is issued
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
    
    def plot(self):
        #t = self.time_vector()
        #dnum = list()
        #for i in range(len(t)):
            #dnum[i] = (t[i].timestamp / 86400) + UTCDateTime(1970, 1, 1).toordinal() 
        snum = self.starttime.timestamp/86400 + UTCDateTime(1970, 1, 1).toordinal() 
        enum = self.endtime.timestamp/86400 + UTCDateTime(1970, 1, 1).toordinal() 
        dnum = np.arange(snum, enum, 1.0/1440)
        plt.figure()
        plt.plot_date(dnum, self.data, '-')
        plt.show()
    
    def time_vector(self):
        t = np.arange(self.starttime, self.endtime, 60)
        return t
