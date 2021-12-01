import os, sys
import datetime as dt
import pprint
from obspy.io.nordic.core import readheader, readwavename, _is_sfile, blanksfile
from obspy.core import read, Stream, UTCDateTime
from obspy.core.event import QuantityError, Pick
from libMVO import correct_nslc

class Sfile:
    'Base class for Sfile parameters'
    
    # Prototypes
    #__init__(self, sfilepath): Construct S-file object from Sfile at path

    def __init__(self, path, use_mvo_parser=False, fast_mode=False):
    
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
        self.error['origintime']=None
        self.error['latitude']=None
        self.error['longitude']=None
        self.error['depth']=None
        self.error['covxy']=None
        self.error['covxz']=None
        self.error['covyz']=None        
        self.focmec = dict()
        self.focmec['strike']=None
        self.focmec['dip']=None
        self.focmec['rake']=None
        self.focmec['agency']=None
        self.focmec['source']=None
        self.focmec['quality']=None     
        self.wavfiles = list()
        self.aeffiles = list()
        self.aefrows = list()
        self.maximum_intensity = None
        self.arrivals = list() 
        self.events = None
        
        
        # Date/time info from sfile path
        self.filetime = spath2datetime(self.path)
        
        if _is_sfile(self.path): 
            if fast_mode:
                self = self.parse_sfile_fast()
            elif use_mvo_parser:
                self = self.parse_sfile()
            else:
                # Here we just parse using ObsPy
                self.events = readheader(self.path)               
                wavfiles = readwavename(self.path) # note this does not include the path
                wavpath = os.path.dirname(self.path).replace('REA','WAV')
                for wavfile in wavfiles:
                    self.wavfiles.append(os.path.join(wavpath, wavfile))

        
    def parse_sfile(self):
        print('Parsing ',self.path)
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
                        if not self.path in _aeffiles:
                            _aeffiles.append(self.path)
                            continue
                            
            if line.find('VOLC MBLYTBH')>-1:
                _aeffiles.append(self.path)
                continue

 
            if len(line) < 80:
                if len(line.strip()) > 0:
                     print("Processing %s: ignoring this line: %s" % (self.path, line, ) )
                continue

            if line[79] == '1':
                
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
                        print('Failed 1-21:\n')
                        print('01233456789'*8)
                        print(line)
                                        
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
                        print('Failed 21-52:\n')
                        print('01233456789'*8)
                        print(line)
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
                        print('Failed 55-64:\n')
                        print('01233456789'*8)
                        print(line)
                if line[63:67].strip():
                    try:
                        self.magnitude.append(float(line[63:67]))
                        #self.magnitude_type.append(magtype_map[line[67]])
                        self.magnitude_type.append(line[67])
                        self.magnitude_agency.append(line[68:71].strip())
                    except:
                        print('Failed 63-72:\n')
                        print('01233456789'*8)
                        print(line)
                if line[71:75].strip():
                    try:
                        self.magnitude.append(float(line[71:75]))
                        #self.magnitude_type.append(magtype_map[line[75]])
                        self.magnitude_type.append(line[75])
                        self.magnitude_agency.append(line[76:79].strip())
                    except:
                        print('Failed 71-80:\n')
                        print('01233456789'*8)
                        print(line)
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
                    print('Failed 18-29:\n')
                    print('01233456789'*8)
                    print(line)
                if aphase == 'AMPL':
                    aamp = line[33:40].strip()
                    aper = line[40:45].strip()
                thisarrival = dict()
                thisarrival['sta'] = asta
                thisarrival['chan'] = achan
                thisarrival['phase'] = aphase
                thisarrival['time'] = atime

                
               
                """
                    thisarrival['time_residual'] = float(line[64:68].strip())
                    thisarrival['weight'] = int(line[68:70].strip())
                    thisarrival['distance_km'] = float(line[72:75].strip())
                    thisarrival['azimuth'] = int(line[77:79].strip())
                """
                if line[64:79].strip():
                    thisarrival['time_residual']=parse_string(line, 64, 68, astype='float', stripstr=True)
                    thisarrival['weight']=parse_string(line, 68, 70, astype='int', stripstr=True)
                    thisarrival['distance_km']=parse_string(line, 72, 75, astype='float', stripstr=True)
                    thisarrival['azimuth']=parse_string(line, 77, 79, astype='int', stripstr=True)
                
                if aphase == 'AMPL':
                    thisarrival['amp'] = aamp
                    thisarrival['per'] = aper
                
                onset='questionable'
                if aphase[0]=='I':
                    onset='impulsive'
                    aphase=aphase[1:]
                elif aphase[0]=='E':
                    onset='emergent'
                    aphase=aphase[1:] 
                    
                if achan:
                    traceID = '.%s..%s' % (asta, achan)
                    Fs = 100
                    if asta[0:2]=='MB' and self.filetime.year < 2005:
                        Fs = 75
                    fixedID = correct_nslc(traceID, Fs, shortperiod=False)
                    #print(traceID,'->',fixedID)                    
                    if 'time_residual' in thisarrival:
                        tres = QuantityError(uncertainty=thisarrival['time_residual'])
                        p = Pick(time = atime, waveform_id = fixedID, onset=onset, phase_hint=aphase, time_errors=tres, backazimuth=180+thisarrival['azimuth'] % 360)
                    else:    
                        p = Pick(time = atime, waveform_id = fixedID, onset=onset, phase_hint=aphase)
                    self.arrivals.append(p)
                
        for _aeffile in _aeffiles:
            self.aeffiles.append(AEFfile(_aeffile))

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
            for pick in self.arrivals:
                str += pick.__str__()
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
    
    def to_dict(self):
        sdict = {}
        sdict['path']=self.path
        sdict['filetime']=self.filetime #.strftime("%Y-%m-%d %H:%M:%S.%f")
        sdict['mainclass']=self.mainclass
        sdict['subclass']=self.subclass
        sdict['wavfile1']=None
        sdict['wavfile2']=None
        
        if len(self.wavfiles)>0:            
            sdict['wavfile1']=self.wavfiles[0].path
           
        if len(self.wavfiles)>1:
            sdict['wavfile2']=self.wavfiles[1].path            

                
        sdict['num_magnitudes']=len(self.magnitude)
        mag,mtype,agency=self.maximum_magnitude()
        sdict['magnitude']=mag
        sdict['magnitude_type']=mtype
        #sdict['magnitude_agency']=agency        

        sdict['num_wavfiles']=len(self.wavfiles)
        sdict['num_aeffiles']=len(self.aeffiles)
 
        
        sdict['located']=False
        if self.longitude:
            sdict['located']=True
        sdict['num_arrivals']=len(self.arrivals)
        #sdict['num_origins']=len(self.arrivals)       
        #sdict['num_events']=len(self.events)
        
        """
        sdict['otime']=self.otime #.strftime("%Y-%m-%d %H:%M:%S.%f")
        sdict['longitude']=self.longitude
        sdict['latitude']=self.latitude
        sdict['depth']=self.depth
        sdict['z_indicator']=self.z_indicator
        sdict['agency']=self.agency
        sdict['id']=self.id
        sdict['rms']=self.rms
        sdict['gap']=self.gap
        sdict['error']=self.error
        """
        if self.error['latitude']:
            sdict['error_exists']=True
        else:
            sdict['error_exists'] = False               
        if self.focmec['strike']:
            sdict['focmec_exists']=True
        else:
            sdict['focmec_exists'] = False
        #sdict['analyst']=self.analyst        

        return sdict
    
    def to_csv(self, csvfile):
        sdict = self.to_dict()
        sdict.to_csv(csvfile)
    
    def to_obspyevent(self):
        thisevent = self.events
        return thisevent
    
    def to_pickle(self):
        import pickle
        picklefile = self.path.replace('REA','PICKLE') + '.pickle'
        with open(picklefile, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)

    def maximum_magnitude(self):
        """ Return the maximum magnitude and corresponding type of magnitude and agency """
        # function mag, mtype, agency = Sfile.maximum_magnitude()
        mag = None
        mtype = None
        agency = None
        for i,magnitude in enumerate(self.magnitude):
            if not mag or magnitude > mag and 'MVO' in self.magnitude_agency: # if not mag needed else will fail when comparing to None
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
        
                 
    def printEvents(self):
        evobj = sfileobj.events
        print(evobj)
        print(evobj.event_descriptions)
        for i, originobj in enumerate(evobj.origins):
            print(i, originobj)
            
    def parse_sfile_fast(self):
        print('Parsing ',self.path)
        fptr = open(self.path,'r') 
        row = fptr.readlines()
        fptr.close()

        # File opened ok, parse each line in a way that depends on the 80th character (line[79])
        for line in row:
            result = line.find("VOLC")
            if line[-1] == '3' or result>-1:
                if line[1:7]== 'ExtMag':
                    continue
                elif line[1:4]=='URL':
                    continue
                elif result:
                    if line[result:result+9]=="VOLC MAIN":
                        self.subclass = line[result+10:result+20].strip()

                            
            if line.find('VOLC MBLYTBH')>-1:
                continue

            if len(line) < 80:
                continue

            if line[79] == '1':
                self.mainclass = line[21:23].strip()
 

            if line[79] == '6':
                wavfiles = line[1:79].split()
                for wavfile in wavfiles:
                    wavfullpath = os.path.join(os.path.dirname(self.path).replace('REA','WAV'), wavfile)
                    self.wavfiles.append(Wavfile(wavfullpath))



            
def spath2datetime(spath):
    # Date/time info from sfile path
    basename = os.path.basename(spath)
    """
    yyyy=int(basename[13:17])
    mm=int(basename[17:19])
    dd=int(basename[0:2])
    hh=int(basename[3:5])
    mi=int(basename[5:7])
    ss=int(basename[8:10])
    return dt.datetime(yyyy,mm,dd,hh,mi,ss)
    """
    fbs = basename.split('.S')
    #print(fbs)
    yyyy = int(fbs[1][0:4])
    mm = int(fbs[1][4:6])
    dd = int(fbs[0][0:2])
    HH = int(fbs[0][3:5])
    MM = int(fbs[0][5:7])
    SS = float(fbs[0][8:10])
    fdt = dt.datetime(yyyy, mm , dd, HH, MM, 0)
    time_change = dt.timedelta(seconds=SS) 
    fdt = fdt + time_change
    return fdt

def filetime2spath(filetime, mainclass='L', db=None, seisan_data=None, fullpath=True):
    # REA/MVOE_/2001/11/19-1951-51L.S200111
    spath = None
    if isinstance(filetime, UTCDateTime) or isinstance(filetime, dt.datetime):
        spath = '%02d-%02d%02d-%02d%s.S%4d%02d' % (filetime.day, filetime.hour, filetime.minute, filetime.second, mainclass, filetime.year, filetime.month)
        if fullpath:
            spath = os.path.join('%4d' % filetime.year, '%02d' % filetime.month, spath)
            if db:
                spath = os.path.join('REA', db, spath)
                if seisan_data:
                    spath = os.path.join(seisan_data, spath)
    return spath

def filetime2wavpath(filetime, mainclass='L', db=None, seisan_data=None, fullpath=True, y2kfix=False, numchans=0):
    # WAV/MVOE_/2002/04/2002-04-27-0321-57S.MVO___014
    wavpath = None
    dbstring = db
    if db=='MVOE_':
        dbstring = 'MVO__'
    if db=='ASNE_':
        dbstring = 'SPN__'
    if len(dbstring)<5:
        dbstring += '_' * (5-len(dbstring))

    if isinstance(filetime, UTCDateTime) or isinstance(filetime, dt.datetime):
        if not y2kfix and filetime.year < 2000:
            wavpath = '%2d%02d-%02d-%02d%02d-%02dS.%s_%03d' % (filetime.year-1900, filetime.month, filetime.day, filetime.hour, filetime.minute, filetime.second, dbstring, numchans)
        else:
            wavpath = '%4d-%02d-%02d-%02d%02d-%02dS.%s_%03d' % (filetime.year, filetime.month, filetime.day, filetime.hour, filetime.minute, filetime.second, dbstring, numchans)
        if fullpath:
            wavpath = os.path.join('%4d' % filetime.year, '%02d' % filetime.month, wavpath)
            if db:
                wavpath = os.path.join('WAV', db, wavpath)
                if seisan_data:
                    wavpath = os.path.join(seisan_data, wavpath)
    return wavpath
            


def parse_string(line, pos0, pos1, astype='float', stripstr=True):
    _s = line[pos0:pos1]
    if stripstr:
        _s = _s.strip()
    if not _s:
        return None
    if astype=='float':           
        return float(_s)
    elif astype=='int':           
        return int(_s)        
    else:
        return _s
    
def read_pickle(picklefile):
    import pickle
    with open('data.pickle', 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        self = pickle.load(f)
    return self



class Wavfile:
    def __init__(self, path=''):
    # class w = Wavfile(path)
        self.path = path.strip()
        self.filetime = None
        self.st = None
        wavbase = os.path.basename(self.path)
        _x = wavbase.split('.')[0]
        _x = _x.split('-')
        if len(_x)==5:
            yyyy = _x[0]
            mm = _x[1]
        elif len(_x)==4:
            yy = _x[0][0:2]
            if yy[0]=='9':
                yyyy = '19' + yy
            elif yy[0]=='2':
                yyyy = '20' + yy
            mm = _x[0][2:4]
        dd = _x[-3]
        HH = _x[-2][0:2]
        MI = _x[-2][2:4]
        SS = _x[-1][0:2]
        #[dd, HH, MI, SS] = _x[-4:0]    
        '''
            mm = _x[5:7]
            dd = _x[8:10]
            HH = _x[11:13]
            MI = _x[13:15]
            SS = _x[16:18]
        elif len(_x)==16:
            yyyy = _x[0:4]
            mm = _x[5:7]
            dd = _x[8:10]
            HH = _x[11:13]
            MI = _x[13:15]
            SS = _x[16:18] 
        '''
        #print(yyyy,mm,dd,HH,MI,SS)    
        self.filetime = dt.datetime(int(yyyy), int(mm), int(dd), hour=int(HH), minute=int(MI), second=int(SS))

    
    def register(self, evtype, userid, overwrite=False, evtime=None):
        if not evtime:
            evtime = self.filetime
        sfilepath = blanksfile(wavfile, evtype, userid, overwrite=False, evtime=evtime)
        return sfilepath
    
    def find_sfile(self, mainclass='L'):
        """ incomplete. idea is to find the sfile, which may not have the same time """
        was_found = False
        """
        spath = os.path.dirname(self.path).replace('WAV','REA',1)
        yyyy = self.filetime.year
        mm = self.filetime.month
        dd = self.filetime.day
        HH = self.filetime.hour
        MI = self.filetime.minute
        SS = self.filetime.second
        sbase = dd + "-" + HH + MI + "-" + SS + "%s.S" % mainclass + yyyy + mm
        sfile = os.path.join(spath, sbase)
        """
        db = os.path.dirname(self.path).split('/')[-3]
        seisan_data = self.path.split('/WAV')[0]
        print(seisan_data, db)
        sfile = filetime2spath(self.filetime, mainclass=mainclass, db=db, seisan_data=seisan_data, fullpath=True)
        if os.path.exists(sfile):
            was_found = True       
        return sfile, was_found
               
    def read(self):
        if os.path.exists(self.path):
            try:
                self.st = read(self.path)
            except:
                self.st = None
                print('failed to read')
                return False
            else:
                return True
        else:
            print('failed to find')
            return False

    def plot(self, equal_scale=False): # must read before plotting
        if self.st:
            self.st.plot(equal_scale=equal_scale)
            return True
        else:
            return False            
            
    def mulplt(self): # must read before plotting
        if self.st:
            libseisGT.mulplt(self.st)
            return True
        else:
            return False


    def __str__(self):
        #str = "wavfile: %s %r " % (self.path, os.path.exists(self.path))
        #str = "wavfile: %s %d " % (self.path, os.path.exists(self.path))
        str = "wavfile: " + (self.path)
        if self.st:
            str += "\n\t" + self.st.__str__()
        return str


class AEFfile:

    def __init__(self, path=None):
        self.aefrows = list() # each element is an AEFrow dict
        self.trigger_window =  None
        self.average_window = None
        self.path = path.strip()
        if not os.path.exists(path):
            print("%s does not exist" % path)
            return
        try:
            file = open(path,'r') 
            lines = file.readlines()
        except IOError:
            print("Error: The specified file does not exist - %s" % (path))
            raise e

        for line in lines:

            if len(line) < 80:
                continue

            if line[79] == '3' or  (line[79] == ' ' and line[1:5]=="VOLC"):
                #print line[1:10]
                if line[1:5]=="VOLC":
                    if line[1:10]=="VOLC MAIN":
                        pass
                    else:
                        aefrow = AEFfile.parse_aefline(line)
                        if aefrow:
                            self.aefrows.append(aefrow)
                _i = line.find("trigger window")
                if _i > -1:
                    _str1 = line[_i:_i+24]
                    _i2 = _str1.find('=') + _i
                    _i3 = _str1.find('s') + _i
                    trigwinstr = line[_i2+1:_i3].strip()
                    self.trigger_window = float(trigwinstr)
                _i = line.find("average window")
                if _i > -1:
                    _str1 = line[_i:_i+24]
                    _i2 = _str1.find('=') + _i
                    _i3 = _str1.find('s') + _i
                    self.average_window = float(line[_i2+1:_i3])
        return None

    def __str__(self):
        #str = "aeffile: %s %r " % (self.path, os.path.exists(self.path))
        str = "aeffile: %s " % self.path
        #str = "\n\texists: " + str(os.path.exists(self.path))
        if self.trigger_window:
            str += "\n\ttrigger window: %f" % self.trigger_window
        if self.average_window:
            str += "\n\taverage window: %f" % self.average_window
        for aefrow in self.aefrows:
            str += aefrow.__str__()
        return str
            
    def cat(self):
        print(cat(self.path))
        return None

    def parse_aefline(line):
        aefrow = {'station':None, 'channel':None, 'amplitude':None, 'energy':None, 'ssam':None, 'maxf':None}
        a_idx = line[15:22].find('A')+15 # where amplitude info starts
        #try:
        station = line[6:10].strip()
        channel = line[11:14].strip()
        
        aefrow['id'] = '.%s..%s' % (station, channel)
        if station[0:2]=='MB':            
            aefrow['fixed_id'] = correct_nslc(aefrow['id'], 100.0, shortperiod=False)
        else:
            aefrow['fixed_id'] = correct_nslc(aefrow['id'], 100.0, shortperiod=True)
        aefrow['amplitude'] = float(line[a_idx+1:a_idx+9].strip())
        aefrow['energy'] = float(line[a_idx+11:a_idx+19].strip())
        aefrow['ssam'] = AEFfile.parse_F(line,aefrow['energy'],a_idx+21)
        if a_idx < 20:
            aefrow['maxf'] = float(line[73:79].strip())
        #except:
        #    print(aefrow)
        return aefrow
    
    def parse_F(line, energy, startindex):
        F = {"frequency_bands":[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 30.0], 
             "percentages": [], "energies": []}
        while startindex < 79 and len(F["percentages"])<12:
            ssamstr = line[startindex:startindex+3].strip()
            if not "." in ssamstr:
                thisval = int(ssamstr)
                F["percentages"].append(thisval)
                F["energies"].append(thisval/100.0 * energy)
            startindex += 3
        return F
