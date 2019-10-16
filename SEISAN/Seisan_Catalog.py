#!/Users/thompsong/miniconda3/bin/python
import os
import datetime
import pprint
import numpy as np
import sys, glob, obspy.core
#SEISAN_DATA_ROOT = ""

def event_list(startdate,enddate):
    event_list=[]
    reapath = os.path.join(SEISAN_DATA_ROOT, 'REA', DB)
    years=list(range(startdate.year,enddate.year+1))
    for year in years:
        #print year
        if year==enddate.year and year==startdate.year:
            months=list(range(startdate.month,enddate.month+1))
        elif year==startdate.year:
            months=list(range(startdate.month,13))
        elif year==enddate.year:
            months=list(range(1,enddate.month+1))
        else:
            months=list(range(1,13))
        for month in months:
            #print month
            dir=os.path.join(reapath, "%04d" % year, "%02d" % month)
            #print dir
            flist=glob(os.path.join(dir,"*L.S*"))
            event_list.extend(flist)
    return event_list 

def cat(path):
    # List the file contents here
    str = ''
    try:
        fptr = open(path,'r')
        row = fptr.readlines()
        for line in row:
            str += line 
    except IOError:
        pass
    return str

def generate_monthly_csv(mm):
    sfileslist = sorted(glob.glob(os.path.join(mm, "*")))
    outfile = mm[-13:-8] + 'catalog' + mm[-7:-3] + mm[-2:] + '.csv'
    print("...Creating %s" % outfile)
    fptr = open(outfile,'w+')
    fptr.write('datetime, mainclass, subclass, duration, wavfilepath, sampling_rate, npts, traceNum, traceID, sfilepath\n') 
    for thissfile in sfileslist:
        print(thissfile)
        s = Sfile(thissfile)
        #print(s)
        numwavfiles = len(s.wavfiles)
        if not s.subclass:
            s.subclass = "_"
        if numwavfiles > 0:
            for thiswavfile in s.wavfiles:
                if not os.path.exists(thiswavfile.path):
                    continue
                st = obspy.read(thiswavfile.path)
                tracenum = 0
                for tr in st:
                     duration = tr.stats.npts / tr.stats.sampling_rate
                     fptr.write("%s,%s,%s,%8.3f,%s,%3d,%6d,%02d,%s,%s\n" % (s.filetime, s.mainclass, \
                         s.subclass, duration, \
                         thiswavfile.path, \
                         tr.stats.sampling_rate, tr.stats.npts, \
                         tracenum, tr.id, thissfile ))
                     tracenum += 1
    fptr.close()


class Sfile:
    'Base class for Sfile parameters'

    def __init__(self, path):
        mytup = path.split("REA/")
        if len(mytup)>1:
           SEISAN_DATA_ROOT = mytup[0]
           mytup2=mytup[1].split('/')
           DB=mytup2[0]
        else:
           SEISAN_DATA_ROOT = "."
        if SEISAN_DATA_ROOT[0]!='/' and SEISAN_DATA_ROOT[0]!='.':
            SEISAN_DATA_ROOT = './' + SEISAN_DATA_ROOT
        if SEISAN_DATA_ROOT[-1]=="/":
            SEISAN_DATA_ROOT = SEISAN_DATA_ROOT[0:-1]
        #print(SEISAN_DATA_ROOT)
        # Initialize optional variable to NULLs
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
        
        # Initialize temporary variables
        _aeffiles = list() # a list of files containing amplitude-energy-frequency rows

        # Open the Sfile and set the filetime
        self.path = path.strip()
        basename = os.path.basename(self.path)
        yyyy=int(basename[13:17])
        mm=int(basename[17:19])
        dd=int(basename[0:2])
        hh=int(basename[3:5])
        mi=int(basename[5:7])
        ss=int(basename[8:10])
        self.filetime = datetime.datetime(yyyy,mm,dd,hh,mi,ss)
        if not os.path.exists(path):
            print("%s does not exist" % path)
            return
        try:
            fptr = open(path,'r') 
            row = fptr.readlines()
        except IOError:
            print("Error: The specified file does not exist - %s" % (self.path))
            raise e

        # File opened ok, parse each line in a way that depends on the 80th character (line[79])
        for line in row:
            result = line.find("VOLC")
            if line[-1] == '3' or result:
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
                        if not path in _aeffiles:
                            _aeffiles.append(path)

            if len(line) < 80:
                print(line + " - IGNORED")
                continue

            if line[79] == '1':
                if len(line[1:20].strip()) >= 14:
                    self.year = int(line[1:5])
                    self.month = int(line[6:8])
                    day = int(line[8:10])
                    hour = int(line[11:13])
                    minute = int(line[13:15])
                    second = 0.0 # sometimes second is not defined at all, so set to 0.0 and then override with value if it exists
                    dummy = line[15:20].strip()
                    if dummy:
                        second = float(line[15:20].strip())
                    if int(second) == 60:
                        minute += 1
                        second -= 60.0
                    #print("second = %f" % second)
                    if second>60: # sometimes second is like 211 which means 2.11 and not 211 seconds. usually it has a decimal point though, e.g. 1.1
                        second=second/100   
                    
                    self.otime = datetime.datetime(self.year, self.month, day, hour, minute, int(second), 1000 * int(  (second - int(second)) * 1000) )
                self.mainclass = line[21:23].strip()
                if line[23:30]!="       ":
                    self.latitude=float(line[23:30])
                if line[30:38]!="        ":
                    self.longitude=float(line[30:38])
                if line[38:43]!="     ":
                    self.depth=float(line[38:43])
                self.z_indicator=line[43].strip()
                self.agency=line[45:48].strip()
                if line[48:51]!="   ":
                    self.no_sta=int(line[48:51])
                else:
                    self.no_sta=0
                if line[51:55]!="    ":
                    self.rms=float(line[51:55])
                if line[55:59] != '    ':
                    self.magnitude.append(float(line[55:59]))
                    #self.magnitude_type.append(magtype_map[line[59]])
                    self.magnitude_type.append(line[59])
                    self.magnitude_agency.append(line[60:63].strip())
                if line[63:67] != '    ':
                    self.magnitude.append(float(line[63:67]))
                    #self.magnitude_type.append(magtype_map[line[67]])
                    self.magnitude_type.append(line[67])
                    self.magnitude_agency.append(line[68:71].strip())
                if line[71:75] != '    ':
                    self.magnitude.append(float(line[71:75]))
                    #self.magnitude_type.append(magtype_map[line[75]])
                    self.magnitude_type.append(line[75])
                    self.magnitude_agency.append(line[76:79].strip())
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
                    wavfullpath = "%s/WAV/%s/%04d/%02d/%s" % (SEISAN_DATA_ROOT, DB, self.year, self.month, wavfile)
                    self.wavfiles.append(Wavfile(wavfullpath))
                    # for each wavfile there is 0 or 1 AEFfile
                    aeffullpath = wavfullpath.replace('WAV', 'AEF')
                    if os.path.exists(aeffullpath):
                        # This aeffile contains AEF information, so store it in list of aeffiles
                        if not aeffullpath in _aeffiles:
                            _aeffiles.append(aeffullpath)

            # Process Type E line, Hyp error estimates
            if line[79] == 'E':
                self.gap=int(line[5:8])
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
                self.otime=datetime.datetime(_yyyy, _mm, _dd, _hh, _mi, _ss, _ms)
                self.latitude=float(line[23:32].strip())
                self.longitude=float(line[33:43].strip())
                self.depth=float(line[44:52].strip())
                self.rms=float(line[53:59].strip())

            if line[79] == 'I':
                self.last_action=line[8:11]
                self.action_time=line[12:26]
                self.analyst = line[30:33]
                self.id = int(line[60:74])

            if line[79] == ' ' and not result:
                if line.strip() == '':
                    continue
                print(line[1:5])
                print(line)
                asta = line[1:5].strip()
                achan = line[5:8].strip()
                aphase = line[8:16].strip()
                ahour = line[18:20].strip()
                aminute = line[20:22].strip()
                asecond = line[22:28].strip()
                if len(asecond)>0:
                    tmplist = asecond.split('.')
                    asecond = tmplist[0]
                    amillisecond = tmplist[1]
                else:
                    asecond=0
                    amillisecond=0
                atime = datetime.datetime(self.year, self.month, day, int(ahour), int(aminute), int(asecond), 1000 * int( amillisecond) )
                if aphase == 'AMPL':
                    aamp = line[33:40].strip()
                    aper = line[40:45].strip()
                thisarrival = dict()
                thisarrival['sta'] = asta
                thisarrival['chan'] = achan
                thisarrival['phase'] = aphase
                thisarrival['time'] = atime
                if aphase == 'AMPL':
                    thisarrival['amp'] = aamp
                    thisarrival['per'] = aper
                self.arrivals.append(thisarrival)

        for _aeffile in _aeffiles:
            self.aeffiles.append(AEFfile(_aeffile))
        return None

    def __str__(self):
        str = "S-file path: " + self.path
        str += "\n\tStart time: " + self.filetime.strftime("%Y-%m-%d %H:%M:%S.%f")
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
            for wavfile in self.wavfiles:
                str += "\n" + wavfile.__str__()
            for aeffile in self.aeffiles:
                str += "\n" + aeffile.__str__()

        return str

    def to_css(self):
        pass

    def maximum_magnitude(self):
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
        print(cat(self.path))
        return None

class Wavfile:
    def __init__(self, path=None):
        self.path = path.strip()
        self.st = None
        return None

    def plot(self):
        if self.st:
            self.st.plot()
        else:
            self.read()
            self.st.plot()
        return True

    def mulplt(self):
        if self.st:
            self.st.mulplt()
        else:
            self.read()
            self.st.mulplt()
        return True

    def read(self):
        #st = None
        self.st = StreamGT()
        if os.path.exists(self.path):
            self.st.read(self.path)
        return True

    def __str__(self):
        #str = "wavfile: %s %r " % (self.path, os.path.exists(self.path))
        #str = "wavfile: %s %d " % (self.path, os.path.exists(self.path))
        str = "wavfile: " + (self.path)
        if self.st:
            str += "\n\t" + self.st.__str__()
        return str

class SSAM:
    def __init__(self, line, energy):
        self.frequency_bands = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 30.0]
        self.percentages = list()
        self.energies = list()
        startindex = 37
        c = 0
        while startindex < 70:
            self.percentages.append(int(line[startindex:startindex+2].strip()))
            self.energies.append(self.percentages[c]/100.0 * energy)
            c += 1
            startindex += 3
        return None
    def __str__(self):
        str = ""
        for i in self.percentages:
            str += "|%d" % i
        str += "|"
        return str

class AEFrow:
    def __init__(self, station, channel, amplitude, energy, ssam, maxf):
        self.station = station
        self.channel = channel
        self.amplitude = amplitude
        self.energy = energy
        self.ssam = ssam
        self.maxf = maxf
        return None

    def __str__(self):
        str = "\n\tstation = %s, channel = %s, amplitude = %e, energy = %e, maxf = %4.2f" \
            % (self.station, self.channel, self.amplitude, self.energy, self.maxf)
        str += ", ssam = " + self.ssam.__str__()
        return str

def parse_aefline(line):
    aefrow = None
    # NEED TO PARSE A LINE 3 HERE FOR ELEMENTS AND THEN CREATE AN AEFROW
    station = line[6:10].strip()
    #print station
    channel = line[11:14].strip()
    #print channel
    amplitude = float(line[16:24])
    energy = float(line[26:34])
    ssam = SSAM(line,energy)
    maxf = float(line[73:78])
    aefrow = AEFrow(station, channel, amplitude, energy, ssam, maxf)
    return aefrow

class AEFfile:

    def __init__(self, path=None):
        self.aefrows = list() # each element is an AEFrow
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
                        aefrow = parse_aefline(line)
                        if aefrow:
                            self.aefrows.append(aefrow)
                _i = line.find("trigger window")
                if _i > -1:
                    _str1 = line[_i:_i+20]
                    _i2 = _str1.find('=') + _i
                    _i3 = _str1.find('s') + _i
                    self.trigger_window = float(line[_i2+1:_i3])
                _i = line.find("average window")
                if _i > -1:
                    _str1 = line[_i:_i+20]
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

if __name__ == "__main__":

    SEISAN_DATA_ROOT = os.path.join("/raid", "data", "seisan")
    DB = "MVOE_"
    sfilelist = [ 
                    '/raid/data/seisan/REA/MVOE_/2002/01/09-1204-48L.S200201', 
                    '/raid/data/seisan/REA/MVOE_/1996/11/30-1254-26L.S199611', 
                    '/raid/data/seisan/REA/MVOE_/1996/11/30-1254-26L.S199611', 
                    '/raid/data/seisan/REA/MVOE_/2000/05/16-2339-12L.S200005', 
                    '/raid/data/seisan/REA/MVOE_/1997/10/11-2228-03L.S199710',
                    '/raid/data/seisan/REA/MVOE_/1996/12/24-2103-55R.S199612',
                    '/raid/data/seisan/REA/MVOE_/1996/12/01-0138-00R.S199612' 
                ]

    for path in sfilelist:

        # load an Sfile into a Sfile object
        s1 = Sfile(path)

        # show the file
        s1.cat() # cat(s1.path) would also work

        # print the object
        print(s1)
    
        # export to ObsPy event object
        #e = s1.to_css()
        #print e

