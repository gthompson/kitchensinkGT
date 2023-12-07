#! /usr/bin/env python
"""
CatalogGT is an extension module for ObsPy Catalog objects
Glenn Thompson, 2014/03/06
"""
import types
import scipy
from CSS_Catalog import CSS_Catalog
import pylab as plt
#from obspy import UTCDateTime
from obspy.core.event import Catalog
#from obspy.core.event import Event, Origin, Magnitude
from obspy.core.util.geodetics import FlinnEngdahl
#import cPickle as pickle
import pickle

def readCSS(dbname):
    """ Call this like:
            import Catalog, CatalogGT
            cat = readCSS(dbname)
    """
    print(CSS_Catalog)
    c = CSS_Catalog()
    c.read(dbname)
    cat = c.to_obspy()
    return cat

def plotlonlat(lon, lat, FontSize=10):
    """ plot latitude versus longitude (Map view) """
    plt.plot(lon, lat, 'o')
    #plt.title('Map view')
    plt.xlabel('Longitude', fontsize=FontSize)
    plt.ylabel('Latitude', fontsize=FontSize)
    #ah = plt.gca()
    #xtl = ah.get_xticklabels()
    #print xtl
    #ah.set_xticklabels(xtl, fontsize=FontSize)

def plotdepthtime(depth, time, FontSize=10):
    """ plot depth against time  """
    plt.plot_date(time, depth,'o')
    #plt.title('Depth-time view')
    plt.xlabel('time', fontsize=FontSize)
    plt.ylabel('depth', fontsize=FontSize)
    plt.gca().invert_yaxis()

def plottimeml(time, ml, FontSize=10):
    """ plot magnitude versus time """
    plt.plot_date(time, ml)
    #plt.title('Mag-time view')
    plt.xlabel('time', fontsize=FontSize)
    plt.ylabel('magnitude', fontsize=FontSize)

def plot_time_ml(ax, time, ml, x_locator, x_formatter, snum, enum):
    """ This is the improved routine from modgiseis """
    ## we do not want to plot Ml = -999.0, the Antelope value for a non-existent Ml - filter them out
    #i = np.where( ml > -3.0)
    #time = time[i]
    #ml = ml[i]

    # plot the data
    ax.plot_date(time, ml, linestyle='o', markerfacecolor='None')
    ax.grid(True)
    ax.xaxis_date()
    plt.setp( ax.get_xticklabels(), rotation=90, horizontalalignment='center', fontsize=7 )
    ax.set_ylabel('Ml', fontsize=8)
    ax.xaxis.set_major_locator(x_locator)
    ax.xaxis.set_major_formatter(x_formatter)
    if snum and enum:
        ax.set_xlim(snum, enum)
    return

def plotdepthlon(depth, lon, FontSize=10):
    """ plot depth vs. longitude """
    plt.plot(lon, depth,'o')
    #plt.title('Depth-lon view')
    plt.xlabel('longitude', fontsize=FontSize)
    plt.ylabel('depth', fontsize=FontSize)
    plt.gca().invert_yaxis()

def plotlatdepth(lat, depth, FontSize=10):
    """ plot latitude vs. depth """
    plt.plot(depth, lat,'o')
    #plt.title('Lat-depth view')
    plt.xlabel('depth', fontsize=FontSize)
    plt.ylabel('latitude', fontsize=FontSize)

def plotcounts(time, binsize, pngfile):
    """ bin the data based on a bin size of binsize days """
    timemin = scipy.floor(min(time)/86400) * 86400
    timemax = scipy.ceil(max(time)/86400) * 86400
    r = [timemin, timemax]
    nbins = (timemax - timemin) / (binsize * 86400)
    h = scipy.histogram(time, nbins, r)
    ecount = h.__getitem__(0)
    ebin = h.__getitem__(1)
    edt = vector_epoch2datetime(ebin)
    plt.plot_date(edt, ecount, linestyle='steps-mid')

""" THESE ARE THE MODGISEIS BINNING AND COUNTS ROUTINES """

def bin_counts(time, bin_edges_in):
    # count the number of "events" in each bin
    # use this if want to produce counts (per bin) but not actually plot them!
    counts, bin_edges_out = np.histogram(time, bin_edges_in)
    return counts

def plot_counts(ax, time, x_locator, x_formatter, bin_edges_in, snum, enum):
    # compute all data needed
    cumcounts = np.arange(1,np.alen(time)+1)
    if len(bin_edges_in) < 2:
        return
    binsize = bin_edges_in[1]-bin_edges_in[0]
    binsize_str = binsizelabel(binsize)

    # plot
    counts, bin_edges_out, patches = ax.hist(time, bin_edges_in, cumulative=False, histtype='bar', color='black', edgecolor=None)
    ax.grid(True)
    ax.xaxis_date()
    plt.setp( ax.get_xticklabels(), rotation=90, horizontalalignment='center', fontsize=7 )
    ax.set_ylabel("# Earthquakes\n%s" % binsize_str, fontsize=8)
    ax.xaxis.set_major_locator(x_locator)
    ax.xaxis.set_major_formatter(x_formatter)
    if snum and enum:
        ax.set_xlim(snum, enum)
    ax2 = ax.twinx()
    p2, = ax2.plot(time,cumcounts,'g', lw=2.5)
    ax2.yaxis.get_label().set_color(p2.get_color())
    ytl_obj = plt.getp(ax2, 'yticklabels')  # get the properties for yticklabels
    #plt.getp(ytl_obj)                       # print out a list of properties
    plt.setp(ytl_obj, color="g")            # set the color of yticks to red
    plt.setp(plt.getp(ax2, 'yticklabels'), color='g') #xticklabels: same
    ax2.set_ylabel("Cumulative\n# Earthquakes", fontsize=8)
    ax2.xaxis.set_major_locator(x_locator)
    ax2.xaxis.set_major_formatter(x_formatter)
    if snum and enum:
        ax2.set_xlim(snum, enum)

""" THESE ARE THE MODGISEIS ENERGY BINNING AND PLOTTING ROUTINES """
def bin_irregular(time, y, bin_edges):
    # bin y against time according to bin_edges (not for binning counts, since they don't have a y value)

    # bin the data as for counts
    counts_per_bin, bin_edges_out = np.histogram(time, bin_edges)
    i_start = 0
    i_end = -1
    binned_y = np.empty(np.alen(counts_per_bin))

    for binnum in range(np.alen(counts_per_bin)):
        i_end += counts_per_bin[binnum]
        if i_start <= i_end:
            binned_y[binnum] = np.sum(y[i_start:i_end+1])
        else:
            binned_y[binnum] = 0
        i_start = i_end + 1
    return binned_y

def ml2energy(ml):
    energy = np.power(10, 1.5 * ml)
    return energy

def energy2ml(energy):
    ml = np.log10(energy)/1.5
    return ml

def plot_energy(ax, time, ml, x_locator, x_formatter, bin_edges, snum, enum):

    # compute all data needed
    energy = ml2energy('ml')
    cumenergy = np.cumsum(energy)
    binned_energy = bin_irregular(time, energy, bin_edges)
    if len(bin_edges) < 2:
        return
    barwidth = bin_edges[1:] - bin_edges[0:-1]
    binsize = bin_edges[1]-bin_edges[0]
    binsize_str = binsizelabel(binsize)

    # plot
    ax.bar(bin_edges[:-1], binned_energy, width=barwidth, color='black', edgecolor=None)

    # re-label the y-axis in terms of equivalent Ml rather than energy
    yticklocs1 = ax.get_yticks()
    ytickvalues1 = np.log10(yticklocs1) / 1.5
    yticklabels1 = list()
    for count in range(len(ytickvalues1)):
        yticklabels1.append("%.2f" % ytickvalues1[count])
    ax.set_yticks(yticklocs1)
    ax.set_yticklabels(yticklabels1)

    ax.grid(True)
    ax.xaxis_date()
    plt.setp( ax.get_xticklabels(), rotation=90, horizontalalignment='center', fontsize=7 )
    ax.set_ylabel("Energy %s\n(unit: Ml)" % binsize_str, fontsize=8)
    ax.xaxis.set_major_locator(x_locator)
    ax.xaxis.set_major_formatter(x_formatter)
    if snum and enum:
        ax.set_xlim(snum, enum)

    # Now add the cumulative energy plot - again with yticklabels as magnitudes
    ax2 = ax.twinx()
    p2, = ax2.plot(time,cumenergy,'g',lw=2.5)

    # use the same ytick locations as for the left-hand axis, but label them in terms of equivalent cumulative magnitude
    yticklocs1 = ax.get_yticks()
    yticklocs2 = (yticklocs1 / max(ax.get_ylim())) * max(ax2.get_ylim() )
    ytickvalues2 = np.log10(yticklocs2) / 1.5
    yticklabels2 = list()
    for count in range(len(ytickvalues2)):
        yticklabels2.append("%.2f" % ytickvalues2[count])
    ax2.set_yticks(yticklocs2)
    ax2.set_yticklabels(yticklabels2)

    ax2.yaxis.get_label().set_color(p2.get_color())
    ytl_obj = plt.getp(ax2, 'yticklabels')  # get the properties for yticklabels
    #plt.getp(ytl_obj)                       # print out a list of properties
    plt.setp(ytl_obj, color="g")            # set the color of yticks to red
    plt.setp(plt.getp(ax2, 'yticklabels'), color='g') #xticklabels: same
    ax2.set_ylabel("Cumulative Energy\n(unit: Ml)",fontsize=8)
    ax2.xaxis.set_major_locator(x_locator)
    ax2.xaxis.set_major_formatter(x_formatter)
    if snum and enum:
        ax2.set_xlim(snum, enum)

""" THE MODGISEIS ROUTINES DEPEND ON COMPUTING bin_edges FIRST """
def autobinsize(daysdiff):
    # Try and keep to around 100 bins or less
    if daysdiff <= 2.0/24:  # less than 2 hours of data, use a binsize of 1 minute
        binsize = 1.0/1440
    elif daysdiff <= 4.0:  # less than 4 days of data, use a binsize of 1 hour
        binsize = 1.0/24
    elif daysdiff <= 100.0:  # less than 100 days of data, use a binsize of 1 day
        binsize = 1.0
    elif daysdiff <= 700.0: # less than 700 days of data, use a binsize of 1 week
        binsize = 7.0
    elif daysdiff <= 365.26 * 23: # less than 23 years of data, use a binsize of (approx) 1 month
        binsize = 365.26/12
    else:
        binsize = 365.26 # otherwise use a binsize of 1 year
    return binsize

def compute_bins(time, snum=None, enum=None, binsize=None):
    # If snum and enum are provided, enum will be end of last bin UNLESS you ask for binsize=365.26, or 365.26/12
    # in which case it will be end of year or month boundary
    # If snum and enum not given, they will end at next boundary - and weeks end on Sat midnight/Sunday 00:00

    # First lets calculate the difference in time between the first and last events
    if (snum==None):
        snum = np.min(time) # time of first event
        enum = np.max(time) # time of last event
    daysdiff = enum - snum

    if (binsize==None):
        binsize = autobinsize(daysdiff)

    # special cases
    if binsize == 365.26/12:
        # because a month isn't exactly 365.26/12 days, this is not going to be the month boundary
        # so let us get the year and the month for snum, but throw away the day, hour, minute, second etc
        sdate = mpl.dates.num2date(snum)
        sdate = datetime.datetime(sdate.year, sdate.month, 1, 0, 0, 0)
        thisyear = sdate.year
        thismonth = sdate.month
        snum = mpl.dates.date2num(sdate)
        bins = list()
        bins.append(snum)
        count = 0
        while bins[count] < enum + binsize:
            count += 1
            thismonth += 1
            if thismonth > 12: # datetime.datetime dies if sdate.month > 12
                thisyear += 1
                thismonth -= 12
            monthdate = datetime.datetime(thisyear, thismonth, 1, 0, 0, 0)
            bins.append(mpl.dates.date2num(monthdate))
        bins = np.array(bins)
        enum = np.max(bins)

    elif binsize == 365.26: # binsize of 1 year
         # because a year isn't exactly 365.26 days, this is not going to be the year boundary
         # so let us get the year for snum, but throw away the month, day, hour, minute, second etc
         sdate = mpl.dates.num2date(snum)
         sdate = datetime.datetime(sdate.year, 1, 1, 0, 0, 0)
         snum = mpl.dates.date2num(sdate)
         bins = list()
         bins.append(snum)
         count = 0
         while bins[count] < enum + binsize:
             count += 1
             yeardate = datetime.datetime(sdate.year + count, 1, 1, 0, 0, 0)
             bins.append(mpl.dates.date2num(yeardate))
         bins = np.array(bins)
         enum = np.max(bins)

    else: # the usual case
        # roundoff the start and end times based on the binsize
        if snum==None and enum==None:
            print("snum and enum undefined - calculating")
            snum = floor(snum, binsize) # start time
            enum = ceil(enum, binsize) # end time
        #bins = np.arange(snum, enum+binsize, binsize)
        bins = np.arange(enum, snum-binsize, -binsize)
        bins = bins[::-1]
    print('snum: %s' % datenum2datestr(snum))
    print('enum: %s' % datenum2datestr(enum))
    return bins, snum, enum

""" FOR LABELLING ONLY - FROM MODGISEIS """
def binsizelabel(binsize):
    binsize_str = ""
    if binsize == 1.0/1440:
        binsize_str = "per minute"
    elif binsize == 1.0/24:
        binsize_str = "per hour"
    elif binsize == 1.0:
        binsize_str = "per day"
    elif binsize == 7.0:
        binsize_str = "per week"
    elif binsize >= 28 and binsize <=31:
        binsize_str = "per month"
    elif binsize >= 365 and binsize <= 366:
        binsize_str = "per year"
    return binsize_str

""" THE FOLLOWING ARE NOT RELATED TO MODGISEIS """


def loadfile(filename):
    with open(filename, 'rb') as input:
        cat = pickle.load(input)
    return cat

def patch(self):
    # Add new methods here that act on catalog pbjects
    """ To add these methods to a Catalog object do:
        patch(cat)
        cat.plotlonlat()
    """
    def volplot(self, pngfile=None):
        """ Replicate Guy's VolPlot routines """
        #lon = list(e.origins[0].longitude for e in self.events)
        #lat = list(e.origins[0].latitude for e in self.events)
        #time = list(e.origins[0].time for e in self.events)
        #depth = list(e.origins[0].depth for e in self.events)
        #mag = list(e.magnitudes[0].mag for e in self.events)
        vd = self.to_vectors()
        FontSize = 12
        plt.close('all')
        plt.figure(1)
        plt.axes([0.1, 0.45, 0.5, 0.5])
        plotlonlat(vd['longitude'], vd['latitude'], FontSize)
        plt.axes([0.1, 0.25, 0.5, 0.15])
        plotdepthlon(vd['depth'], vd['longitude'], FontSize)
        plt.axes([0.1, 0.05, 0.8, 0.15])
        plotdepthtime(vd['depth'], vd['time'], FontSize)
        plt.axes([0.7, 0.45, 0.2, 0.5])
        plotlatdepth(vd['latitude'], vd['depth'], FontSize)
        if pngfile:
            plt.savefig(pngfile)
        else:
            plt.show()
    self.volplot = types.MethodType(volplot, self)

    def to_vectors(self):
        # lat/lon coordinates, magnitudes, dates
        lats = []
        lons = []
        mags = []
        depths = []
        times = []
        for event in self:
            if not event.origins:
                continue

            origin = event.preferred_origin() or event.origins[0]
            lats.append(origin.latitude)
            lons.append(origin.longitude)
            times.append(origin.time)
            depths.append(origin.depth)
            if event.magnitudes:
                magnitude = event.preferred_magnitude() or event.magnitudes[0]
                mag = magnitude.mag
            else:
                mag = -99.99
            mags.append(mag)
        return { 'latitude': lats, 'longitude': lons, 'depth': depths, 'time' : times, 'mag' : mags }
    self.to_vectors = types.MethodType(to_vectors, self)

    """ This is from dbplotevents.xpy """
    def plot_time_counts_energy(self):
        vd = self.to_vectors()
        numevents = len(vd['lat'])
        nummagnitudes = len(vd['mag'])
        if numevents > 0:

            # Let matplotlib automatically decide where to put date (x-axis) tick marks, 
            # and what style of labels to use
            locator = mpl.dates.AutoDateLocator()
            formatter = mpl.dates.AutoDateFormatter(locator)

            # create the figure canvas
            fig1 = plt.figure()
 
            numplots = 3
            if nummagnitudes <= 1:
                numplots = 1
            plotnum = 1

            # add subplot - mag versus time
            if nummagnitudes > 1:
                ax1 = fig1.add_subplot(numplots, 1, plotnum)
                plotnum += 1
                plot_time_ml(ax1, vd['time'], vd['mag'], locator, formatter, snum, enum)

            if numevents > 1:

                # Compute bin_edges based on the first and last event times
                bin_edges, snum, enum = compute_bins(vd['time'], snum, enum)

                # add subplot - counts versus time
                ax2 = fig1.add_subplot(numplots, 1, plotnum)
                plot_counts(ax2, vd['time'], locator, formatter, bin_edges, snum, enum)
                plotnum += 1

                # add subplot - energy versus time
                if nummagnitudes > 1:
                    ax3 = fig1.add_subplot(numplots, 1, plotnum)
                    plotnum += 1
                    plot_energy(ax3, vd['time'], vd['mag'], locator, formatter, bin_edges, snum, enum)

    # Could also add methods here to export to versions of XML and KML used in Alaska
    # and also export to Antelope CSS3.0 again
    # def to_xml(self, target='kml'):
    #     for e in self.events:
    #     if target=='kml':
    #         print '...'

    # Also merge the importer for Seisan

if __name__ == "__main__":
    # Load catalog
    import os
    picklefile = '/tmp/my_custom_events.pk'
    matfile = '/tmp/my_custom_events.mat'
    xmlfile = "/tmp/my_custom_events.xml"
    if os.path.exists(picklefile):
        with open(picklefile, 'rb') as input:
            cat = pickle.load(input)
    elif os.path.exists(matfile):
        from scipy.io import loadmat
        mydict = loadmat(matfile)
        cat = mydict['cat']
    elif os.path.exists(xmlfile):
        from obspy.core.event import readEvents
        cat = readEvents(xmlfile)
    else:
        # read in the Antelope demo database into a CSS_Catalog object
        # and convert to ObsPy Catalog object
    print("Reading in Antelope demo database")
        dbname = '/opt/antelope/data/db/demo/demo'
        cat = readCSS(dbname)

    # Save catalog
    if not os.path.exists(picklefile):
        with open(picklefile, 'wb') as output:
            pickle.dump(cat, output, pickle.HIGHEST_PROTOCOL)
    if not os.path.exists(matfile):
        mydict = {'cat': cat}
        try:
            from scipy.io import savemat
            savemat(matfile, mydict)
        except:
            # TypeError: Could not convert None (type <type 'NoneType'>) to array
            pass
    if not os.path.exists(xmlfile):
        # write the ObsPy Catalog object to QuakeML
        cat.write(xmlfile, format="QUAKEML")
  
    # Plot the events
    #cat.plot(projection="ortho", resolution="l")
    ##cat.plot(projection="local", resolution="i")

    # used some of the plotters defined here
    patch(cat)
    #cat.volplot()
    cat.plotcounts()
    plt.show()


