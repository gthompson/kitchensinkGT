#! /usr/bin/env python
""" sam class, Glenn Thompson 2010/05/01 """
""" Goal is to be able to have an M*N matrix of sam objects - and then plot them accordingly 
So need to init an M*N matrix of sam objects 
Then load an M*N matrix
Filter them
Plot them
"""

class sam(object):
    
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
        print('%s' % self.name,' %12d' % self.estart,' %12d' % self.eend,' %s' % self.sta, ' %s' % self.chan, ' %s' % self.filepattern)
            
    def import_samfile(self):
    """
    open file
    work out number of samples to get
    position pointer
    read samples
    close file
    return data
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

 