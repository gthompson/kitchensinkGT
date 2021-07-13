import os
import datetime as dt
import pprint
from obspy.io.nordic.core import readheader, readwavename, _is_sfile
from Wavfile import Wavfile

class Sfile:
    'Base class for Sfile parameters'
    
    # Prototypes
    #__init__(self, sfilepath): Construct S-file object from Sfile at path

    def __init__(self, path):
    
        self.path = path.strip()
        if not os.path.exists(self.path):
            return
        
        # Initialize optional variable to NULLs
        self.filetime = None # get year as self.filetime.year but might need to convert int to str
        self.otime = None # origin time, a datetime
        self.mainclass = None
        self.latitude = None
        self.longitude = None
        self.depth = None
        self.z_indicator = None
        self.agency = None
        self.no_sta = None
        self.rms = None
        self.magnitude = list()
        self.magnitude_type = list()
        self.magnitude_agency = list()
        self.last_action = None
        self.action_time = None
        self.analyst = None
        self.id = None
        self.subclass = None
        self.url = None
        self.gap = None
        self.error = dict()
        self.focmec = dict()
        self.wavfiles = list()
        self.aeffiles = list()
        self.aefrows = list()
        self.maximum_intensity = None
        self.arrivals = list() 
        self.events = None
        
        
        # Date/time info from sfile path
        basename = os.path.basename(self.path)
        yyyy=int(basename[13:17])
        mm=int(basename[17:19])
        dd=int(basename[0:2])
        hh=int(basename[3:5])
        mi=int(basename[5:7])
        ss=int(basename[8:10])
        self.filetime = dt.datetime(yyyy,mm,dd,hh,mi,ss)
        #self.otime = self.filetime
        
        if _is_sfile(self.path):
             self.events = readheader(self.path)
             #self.wavfiles = readwavename(self.path) # note this does not include the path
             
    def read_event(self):

        fptr = open(self.path,'r') 
        row = fptr.readlines()
        fptr.close()
        
        # Initialize temporary variables
        _aeffiles = list() # a list of files containing amplitude-energy-frequency rows

        # File opened ok, parse each line in a way that depends on the 80th character (line[79])
        for line in row:
            result = line.find("VOLC")
            if line[-1] == '3' or result>-1:
                if line[1:7]== 'ExtMag':
                    self.magnitude.append(float(line[8:12]))
                    #self.magnitude_type.append(magtype_map[line[12]])
                    self.magnitude_type.append(line[12])
                    self.magnitude_agency.append(line[13:16])
                elif line[1:4]=='URL':
                    self.url=line[5:78].rstrip()
                elif result:
                    if line[result:result+9]=="VOLC MAIN":
                        self.subclass = line[result+10:result+20].strip()
                    else:
                        # We think this Sfile contains AEF information, so store it in list of aeffiles
                        if not self.sfilepath in _aeffiles:
                            _aeffiles.append(self.sfilepath)
                            continue
                            
            if line.find('VOLC MBLYTBH')>-1:
                _aeffiles.append(self.sfilepath)
                continue

 
            if len(line) < 80:
                if len(line.strip()) > 0:
                    print("Processing %s: ignoring this line: %s" % (path, line, ) )
                continue

            if line[79] == '1':
                print('Processing: %s' % line)
                if len(line[1:20].strip()) >= 14:
                    try:
                        oyear = int(line[1:5])
                        omonth = int(line[6:8])
                        oday = int(line[8:10])
                        ohour = int(line[11:13])
                        ominute = int(line[13:15])
                        osecond = 0.0 # sometimes second is not defined at all, so set to 0.0 and then override with value if it exists
                        dummy = line[15:20].strip()
                        if dummy:
                            osecond = float(line[15:20].strip())
                        if int(osecond) == 60:
                            ominute += 1
                            osecond -= 60.0
                        #print("second = %f" % second)
                        if osecond>60: # sometimes second is like 211 which means 2.11 and not 211 seconds. usually it has a decimal point though, e.g. 1.1
                            osecond=osecond/100   
                            self.otime = dt.datetime(oyear, omonth, oday, ohour, ominute, int(osecond), 1000 * int(  (osecond - int(osecond)) * 1000) )
                    except:
                        print('Failed:\n', line)
                        pass
                                        
                self.mainclass = line[21:23].strip()
                if line[23:30].strip():
                    self.latitude=float(line[23:30])
                if line[30:38].strip():
                    self.longitude=float(line[30:38])
                if line[38:43].strip():
                    self.depth=float(line[38:43])
                self.z_indicator=line[43].strip()
                self.agency=line[45:48].strip()
                nostastr = line[49:51].strip()
                if nostastr:
                    try:
                        self.no_sta=int(nostastr)
                    except:
                        print('Failed:\n', line)
                        pass
                else:
                    self.no_sta=0
                if line[51:55].strip():
                    self.rms=float(line[51:55])
                if line[55:59].strip():
                    try:
                        self.magnitude.append(float(line[55:59]))
                        #self.magnitude_type.append(magtype_map[line[59]])
                        self.magnitude_type.append(line[59])
                        self.magnitude_agency.append(line[60:63].strip())
                    except:
                        print('Failed:\n', line)
                        pass
                if line[63:67].strip():
                    try:
                        self.magnitude.append(float(line[63:67]))
                        #self.magnitude_type.append(magtype_map[line[67]])
                        self.magnitude_type.append(line[67])
                        self.magnitude_agency.append(line[68:71].strip())
                    except:
                        print('Failed:\n', line)
                        pass
                if line[71:75].strip():
                    try:
                        self.magnitude.append(float(line[71:75]))
                        #self.magnitude_type.append(magtype_map[line[75]])
                        self.magnitude_type.append(line[75])
                        self.magnitude_agency.append(line[76:79].strip())
                    except:
                        print('Failed:\n', line)
                        pass
#                else: #If there is an additional line1 process the additional magnitudes
#                    if line[55:59] != '    ':
#                        self.magnitude.append(float(line[55:59]))
#                        #self.magnitude_type.append(magtype_map[line[59]])
#                        self.magnitude_type.append(line[59])
#                        self.magnitude_agency.append(line[60:63].strip())
#                    if line[63:67] != '    ':
#                        self.magnitude.append(float(line[63:67]))
#                        #self.magnitude_type.append(magtype_map[line[67]])
#                        self.magnitude_type.append(line[67])
#                        self.magnitude_agency.append(line[68:71].strip())
#                    if line[71:75] != '    ':
#                        self.magnitude.append(float(line[71:75]))
#                        #self.magnitude_type.append(magtype_map[line[75]])
#                        self.magnitude_type.append(line[75])
#                        self.magnitude_agency.append(line[76:79].strip())

            # Process Type 2 line, Macroseismic Intensity Information
            if line[79] == '2':
               self.maximum_intensity=int(line[27:29]) 

            if line[79] == '6':
                #print("got a WAV file")
                wavfiles = line[1:79].split()
                #wavfile = line[1:36]
                for wavfile in wavfiles:
                    wavfullpath = os.path.join(os.path.dirname(self.path).replace('REA','WAV'), wavfile)
                    self.wavfiles.append(Wavfile(wavfullpath))
                    # for each wavfile there is 0 or 1 AEFfile
                    aeffullpath = wavfullpath.replace('WAV', 'AEF')
                    if os.path.exists(aeffullpath):
                        # This aeffile contains AEF information, so store it in list of aeffiles
                        if not aeffullpath in _aeffiles:
                            _aeffiles.append(aeffullpath)

            # Process Type E line, Hyp error estimates
            if line[79] == 'E':
                gapstr = line[5:8].strip()
                if len(gapstr)>0:
                    self.gap=int(gapstr)
                else:
                    self.gap = -1
                self.error['origintime']=float(line[14:20])
                self.error['latitude']=float(line[24:30].strip())
                self.error['longitude']=float(line[32:38].strip())
                self.error['depth']=float(line[38:43].strip())
                self.error['covxy']=float(line[43:55].strip())
                self.error['covxz']=float(line[55:67].strip())
                self.error['covyz']=float(line[67:79].strip())

            # Process Type F line, Fault plane solution
            # Format has changed need to fix AAH - 2011-06-23
            if line[79] == 'F' and 'dip' not in self.focmec:
                self.focmec['strike']=float(line[0:10])
                self.focmec['dip']=float(line[10:20])
                self.focmec['rake']=float(line[20:30])
                #self.focmec['bad_pols']=int(line[60:66])
                self.focmec['agency']=line[66:69]
                self.focmec['source']=line[70:77]
                self.focmec['quality']=line[77]

            # Process Type H line, High accuracy line
            # This replaces some origin parameters with more accurate ones
            if line[79] == 'H':
                _osec=float(line[16:22])
                _yyyy=int(line[1:5])
                _mm=int(line[6:8])
                _dd=int(line[8:10])
                _hh=int(line[11:13])
                _mi=int(line[13:15])
                _ss=int(_osec)
                _ms=int((_osec-int(_osec))*1.e6)
                self.otime=dt.datetime(_yyyy, _mm, _dd, _hh, _mi, _ss, _ms)
                self.latitude=float(line[23:32].strip())
                self.longitude=float(line[33:43].strip())
                self.depth=float(line[44:52].strip())
                self.rms=float(line[53:59].strip())

            if line[79] == 'I':
                self.last_action=line[8:11]
                self.action_time=line[12:26]
                self.analyst = line[30:33]
                self.id = int(line[60:74])

            #if line[79] == ' ' and not result:
            if line[79] == ' ' and line[1] == 'M':
                if line.strip() == '':
                    continue
                asta = line[1:5].strip()
                achan = line[5:8].strip()
                aphase = line[8:16].strip()
                ahour = int(line[18:20].strip())
                aminute = int(line[20:22].strip())
                asecond = line[22:28].strip()
                if len(asecond)>0:
                    tmplist = asecond.split('.')
                    asecond = int(tmplist[0])
                    amillisecond = int(tmplist[1])
                else:
                    asecond=0
                    amillisecond=0
                if asecond >= 60:
                    aminute += 1
                    asecond -= 60
                    if aminute >= 60:
                        aminute -= 60
                        ahour += 1
                if ahour > 23:
                    ahour -= 24
                try:
                    atime = self.filetime
                    if self.otime:
                        atime = self.otime
                    ayear = atime.year
                    amonth = atime.month
                    aday = atime.day
                    atime = dt.datetime(ayear, amonth, aday, ahour, aminute, asecond, 1000 * amillisecond)
                except:
                    pass
                if aphase == 'AMPL':
                    aamp = line[33:40].strip()
                    aper = line[40:45].strip()
                thisarrival = dict()
                thisarrival['sta'] = asta
                thisarrival['chan'] = achan
                thisarrival['phase'] = aphase
                thisarrival['time'] = atime
                try:
                    thisarrival['tres'] = float(line[64:68].strip())
                    thisarrival['w'] = int(line[68:70].strip())
                    thisarrival['dis'] = float(line[72:75].strip())
                    thisarrival['caz'] = int(line[77:79].strip())
                except:
                    print('Failed:\n', line)
                    pass
                if aphase == 'AMPL':
                    thisarrival['amp'] = aamp
                    thisarrival['per'] = aper
                self.arrivals.append(thisarrival)
        # Add code ehre to use new AEFfile class to read AEF files
        #for _aeffile in _aeffiles:
            #self.aeffiles.append(AEFfile(_aeffile))


    def __str__(self):
        """ summary string of an Sfile object, for print command """
        str = "S-file path: " + self.path
        str += "\n\tFile time: " + self.filetime.strftime("%Y-%m-%d %H:%M:%S.%f")
        if self.otime:
            str += "\n\tOrigin time: " + self.otime.strftime("%Y-%m-%d %H:%M:%S.%f")
        if self.mainclass:
            str += "\n\tMainclass: " + self.mainclass
        if self.subclass:
            str += "\n\tSubclass: " + self.subclass
        if self.longitude:
            str += "\n\tLongitude: %f " %  self.longitude
        if self.latitude:
            str += "\n\tLatitude: %f" % self.latitude
        if self.depth:
            str += "\n\tDepth: %f" % self.depth
        if self.z_indicator:
            str += "\n\tZ_indicator: " + self.z_indicator
        if self.agency:
            str += "\n\tAgency: " + self.agency
        if self.id:
            str += "\n\tId: %d " % self.id
        if self.rms:
            str += "\n\trms: %f" % self.rms
        for count in range(len(self.magnitude)):
            str += "\n\tMagnitude: %.1f, %s, %s" % (self.magnitude[count], self.magnitude_type[count], self.magnitude_agency[count])
        if self.last_action:
            str += "\n\tLast action: " + self.last_action
        if self.action_time:
            str += "\n\tLast action time: " + self.action_time
        if self.analyst:
            str += "\n\tAnalyst: " + self.analyst
        if self.url:
            str += "\n\tURL: " + self.url
        if self.gap:
            str += "\n\tGap: %d" % self.gap
        if self.error:
            str += "\n\tError: " + pprint.pformat(self.error)
        if self.focmec:
            str += "\n\tFocmec: " + pprint.pformat(self.focmec)
        str += "\n\tNumber of arrivals: %d" % len(self.arrivals)
        if len(self.arrivals)>0:
            str += "\n\tArrivals:" 
            for arrival in self.arrivals:
                str += "\n\t\t" + "%s.%s %s " % (arrival['sta'], arrival['chan'], arrival['phase']) + arrival['time'].strftime("%Y-%m-%d %H:%M:%S.%f")
        if len(self.wavfiles)>0:
            str += "\n\tWAV files:" 
            for wavfile in self.wavfiles:
                str += "\n\t\t" + wavfile.__str__()
        if len(self.aeffiles)>0:
            str += "\n\tAEF files:" 
            for aeffile in self.aeffiles:
                str += "\n\t\t" + aeffile.__str__()
        if self.events:
            str += "\n\tEvents:"
            for event in self.events:
                str += "\n\t\t" + event.__str__()

        return str

    def to_css(self):
        """ Write CSS3.0 database event from an Sfile
        This is potentially less lossy than via an ObsPy Catalog object first """
        pass
    
    def to_obspyevent(self):
        thisevent = self.events
        return thisevent

    def maximum_magnitude(self):
        """ Return the maximum magnitude and corresponding type of magnitude and agency """
        # function mag, mtype, agency = Sfile.maximum_magnitude()
        mag = 0.0
        mtype = "NA"
        agency = "NA"
        for i,magnitude in enumerate(self.magnitude):
            if magnitude > mag:
                mag = magnitude
                mtype = self.magnitude_type[i]
                agency = self.magnitude_agency[i]
        return mag,mtype,agency 

    def cat(self):
        """ Show contents of an Sfile """
        fptr = open(self.path,'r')
        contents = fptr.read()
        fptr.close()
        print(contents)
