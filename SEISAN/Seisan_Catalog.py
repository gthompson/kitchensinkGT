#!/usr/bin/env python
import os
import pprint
import numpy as np
import sys, glob, obspy.core
import pandas as pd
import csv
import trace_quality_control as qc
import datetime as dt
import obspy
SEISAN_DATA = os.environ['SEISAN_DATA']
if not SEISAN_DATA:
    SEISAN_DATA = "./seismo"

def create_event_list(startdate,enddate):
# function event_list = create_event_list(startdate, enddate)
    event_list=[]
    reapath = os.path.join(SEISAN_DATA, 'REA', DB)
    years=list(range(startdate.year,enddate.year+1))
    for year in years:
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
            flist=glob(os.path.join(dir,"*L.S*"))
            event_list.extend(flist)
    return event_list 

def cat(path):
# function str = cat(file)
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

def generate_monthly_csv(mm, flag_sfiles_only=False):
# function generate_monthly_csv(monthdirs, flag_sfiles_only)
    print("Reading %s" % mm)
    sfileslist = sorted(glob.glob(os.path.join(mm, "*")))
    if flag_sfiles_only:
        outfile = mm[-13:-8] + 'sfiles' + mm[-7:-3] + mm[-2:] + '.csv'
        hypofile = mm[-13:-8] + 'hypo' + mm[-7:-3] + mm[-2:] + '.csv'
    else:
        outfile = mm[-13:-8] + 'catalog' + mm[-7:-3] + mm[-2:] + '.csv'
    print("...Creating %s" % outfile)
    if os.path.exists(outfile):
        print("CSV file %s already exists. Not overwriting. If you want to overwrite then please delete from the command line" % outfile)
        return
    fptr = open(outfile,'w')
    if flag_sfiles_only:
        fptr.write('datetime,mainclass,subclass,sfilepath,analyst,wavfile1,wavfile2,numarrivals,nummagnitudes,numaefrows,duration\n') 
    else:
        fptr.write('datetime,mainclass,subclass,duration,wavfilepath,sampling_rate,npts,traceNum,traceID,sfilepath,analyst,fixedID,quality_factor,snr,highval,lowval\n') 
    for thissfile in sfileslist:
        s = Sfile(thissfile)
        if not np.isnan(s.latitude) and flag_sfiles_only:
            if not s.longitude:
                s.longitude = float('nan')
            if not s.depth:
                s.depth = float('nan')
            # write to hypo file
            hfptr = open(hypofile,'a+')
            if not os.path.exists(hypofile):
                hfptr.write('datetime,mainclass,subclass,latitude,longitude,depth,magnitude,magnitude_type,sfilepath\n') 
            if len(s.magnitude)>0:
                magnitude = s.magnitude[0]
                magnitude_type = s.magnitude_type[0]
            else:
                magnitude = float('nan')
                magnitude_type = ""
            try:
                hfptr.write("%s,%s,%s,%9.4f,%9.4f,%6.1f,%4.1f,%s,%s\n" % (s.filetime, s.mainclass, \
                    s.subclass, \
                    s.latitude, s.longitude, s.depth, magnitude, magnitude_type, \
                    thissfile ))
            except:
                print(s.latitude)
                print(s.longitude)
                print(s.depth)
                print(magnitude)
                barf
            hfptr.close()
        numwavfiles = len(s.wavfiles)
        tmp = thissfile.split('REA')
        shortsfile = 'REA' + tmp[1] 
        if not s.subclass:
            s.subclass = "_"
        if flag_sfiles_only:
            wavfile1 = "-"
            wavfile2 = "-"
            if numwavfiles>0:
                tmp = s.wavfiles[0].path.split('WAV')
                wavfile1 = 'WAV' + tmp[1]
            if numwavfiles>1:
                tmp = s.wavfiles[1].path.split('WAV')
                wavfile2 = 'WAV' + tmp[1]
            if len(s.aeffiles)>0:
                duration = s.aeffiles[0].trigger_window
            else:
                duration = 0.0
            if not duration:
                duration = 0.0

            try:
                fptr.write("%s,%s,%s,%s,%s,%s,%s,%2d,%1d,%2d,%5.1f\n" % (s.filetime, s.mainclass, \
                    s.subclass, \
                    thissfile, s.analyst, wavfile1, wavfile2, \
                    len(s.arrivals), len(s.magnitude), len(s.aefrows), duration  ))
            except:
                print(len(s.arrivals))
                print(len(s.magnitude))
                print(len(s.aefrows))
                print(duration)
                barf

        else:
            for thiswavfile in s.wavfiles:
                if not os.path.exists(thiswavfile.path):
                    continue
                try:
                    st = obspy.read(thiswavfile.path)
                except:
                    print("Processing %s: Cannot read wavfile %s" % (thissfile, thiswavfile.path,) )
                    continue
                tracenum = 0
                tmp = thiswavfile.path.split('WAV')
                wavfile1 = 'WAV' + tmp[1]
                for tr in st:
                     tr2, quality_factor, snr = qc.compute_metrics(tr)
                     duration = tr.stats.npts / tr.stats.sampling_rate
                     fptr.write("%s,%s,%s,%8.3f,%s,%3d,%6d,%02d,%s,%s,%s,%s,%d,%.2f,%.2f,%.2f\n" % (s.filetime, s.mainclass, \
                         s.subclass, duration, \
                         wavfile1, \
                         tr.stats.sampling_rate, tr.stats.npts, \
                         tracenum, tr.id, shortsfile, s.analyst, tr2.id, quality_factor, snr[0], snr[1], snr[2] ))
                     tracenum += 1
            
    fptr.close()


class Sfile:
    'Base class for Sfile parameters'

    def __init__(self, path):
    # class s=Sfile(file)
        mytup = path.split("REA/")
        global SEISAN_DATA
        if len(mytup)>1:
           SEISAN_DATA = mytup[0]
           mytup2=mytup[1].split('/')
           DB=mytup2[0]
        else:
           SEISAN_DATA = "."
        if SEISAN_DATA[0]!='/' and SEISAN_DATA[0]!='.':
            SEISAN_DATA = './' + SEISAN_DATA
        if SEISAN_DATA[-1]=="/":
            SEISAN_DATA = SEISAN_DATA[0:-1]
        #print(SEISAN_DATA)
        # Initialize optional variable to NULLs
        self.otime = None # origin time, a datetime
        self.mainclass = None
        self.latitude = float('nan')
        self.longitude = float('nan')
        self.depth = float('nan')
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
        self.filetime = dt.datetime(yyyy,mm,dd,hh,mi,ss)
        self.otime = self.filetime
        self.year = yyyy
        self.month = mm
        self.day = dd
        if not os.path.exists(path):
            print("%s does not exist" % path)
            return
        try:
            fptr = open(path,'r') 
            row = fptr.readlines()
            #os.system("cat %s" % path)
        except IOError:
            print("Error: The specified file does not exist - %s" % (self.path))
            raise e

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
                        if not path in _aeffiles:
                            _aeffiles.append(path)
                            continue
            if line.find('VOLC MBLYTBH')>-1:
                _aeffiles.append(path)
                continue

 
            if len(line) < 80:
                if len(line.strip()) > 0:
                    print("Processing %s: ignoring this line: %s" % (path, line, ) )
                continue

            if line[79] == '1':
                if len(line[1:20].strip()) >= 14:
                    try:
                        self.year = int(line[1:5])
                        self.month = int(line[6:8])
                        self.day = int(line[8:10])
                    except:
                        print(line)
                        barf
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
                    
                    self.otime = dt.datetime(self.year, self.month, self.day, hour, minute, int(second), 1000 * int(  (second - int(second)) * 1000) )
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
                        print(line)
                        print(self.otime)
                        print(self.longitude, self.latitude, self.depth)
                        print(nostastr)
                        barf
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
                        pass
                if line[63:67].strip():
                    try:
                        self.magnitude.append(float(line[63:67]))
                        #self.magnitude_type.append(magtype_map[line[67]])
                        self.magnitude_type.append(line[67])
                        self.magnitude_agency.append(line[68:71].strip())
                    except:
                        print("Processing %s" %(path,) )
                        print(line)
                        print(self.magnitude)
                        print(self.magnitude_type)
                        print(self.magnitude_agency)
                        print("There probably was not a second magnitude")
                if line[71:75].strip():
                    try:
                        self.magnitude.append(float(line[71:75]))
                        #self.magnitude_type.append(magtype_map[line[75]])
                        self.magnitude_type.append(line[75])
                        self.magnitude_agency.append(line[76:79].strip())
                    except:
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
                    wavfullpath = "%s/WAV/%s/%04d/%02d/%s" % (SEISAN_DATA, DB, self.year, self.month, wavfile)
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
                    atime = dt.datetime(self.year, self.month, self.day, ahour, aminute, asecond, 1000 * amillisecond)
                except:
                    print(self.year, self.month, self.day, ahour, aminute, asecond)
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
                    pass
                if aphase == 'AMPL':
                    thisarrival['amp'] = aamp
                    thisarrival['per'] = aper
                #print(thisarrival)
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
        if len(self.wavfiles)>0:
            str += "\n\tWAV files:" 
            for wavfile in self.wavfiles:
                str += "\n\t\t" + wavfile.__str__()
        if len(self.aeffiles)>0:
            str += "\n\tAEF files:" 
            for aeffile in self.aeffiles:
                str += "\n\t\t" + aeffile.__str__()

        return str

    def to_css(self):
        pass

    def maximum_magnitude(self):
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
    # function Sfile.cat()
        print(cat(self.path))
        return None

class Wavfile:
    def __init__(self, path=None):
    # class w = Wavfile(file)
        self.path = path.strip()
        self.st = None
        return None

    def plot(self):
    # function Wavfile.plot()
        if self.st:
            self.st.plot()
        else:
            self.read()
            self.st.plot()
        return True

    def mulplt(self):
    # function Wavfile.mulplt()
        if self.st:
            self.st.mulplt()
        else:
            self.read()
            self.st.mulplt()
        return True

    def read(self):
    # function Wavfile.read()
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
    def __init__(self, line, energy, startindex):
    # class SSAM(line, energy, startindex)
        self.frequency_bands = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 30.0]
        self.percentages = list()
        self.energies = list()
        c = 0
        #print(line[startindex:])
        while startindex < 79 and len(self.percentages)<12:
            ssamstr = line[startindex:startindex+3].strip()
            #print(ssamstr) 
            if not "." in ssamstr:
                thisval = int(ssamstr)
                self.percentages.append(thisval)
                self.energies.append(thisval/100.0 * energy)
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
    # class AEFrow(station, channel, amplitude, energy, ssam, maxf)
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
# function AEFrow = parse_aefline(line)
    aefrow = None
    #print(line)
    a_idx = line[15:22].find('A')+15 # where amplitude info starts
    # NEED TO PARSE A LINE 3 HERE FOR ELEMENTS AND THEN CREATE AN AEFROW
    try:
        station = line[6:10].strip()
        channel = line[11:14].strip()
        amplitude = float(line[a_idx+1:a_idx+9].strip())
        energy = float(line[a_idx+11:a_idx+19].strip())
        ssam = SSAM(line,energy,a_idx+21)
        if a_idx < 20:
            maxf = float(line[73:79].strip())
        else:
            maxf = 0.0
    except:
        print(line)
        print(a_idx)
        print(station)
        print(channel)
        print(amplitude)
        print(energy)
        print(ssam)
        print(maxf)
        barf
    aefrow = AEFrow(station, channel, amplitude, energy, ssam, maxf)
    return aefrow

class AEFfile:

    def __init__(self, path=None):
    # class AEFfile(file)
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


def sfilecsv_daycount(list_of_csv_files):
# function sfilecsv_daycount(list_of_csv_files)

    # Combine all the CSV files, and then summarize
    print('\n*************************')
    print('**** Overall summary ****')
    frames = []
    for csvfile in list_of_csv_files:
        df = pd.read_csv(csvfile)
        df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
        if not isinstance(frames, pd.DataFrame):
            frames = df
        else:
            df = df.reset_index(drop=True)
        frames = pd.concat([frames, df], axis=0)

    # how many of each mainclass?
    class_hash = frames.mainclass.value_counts()
    for ck in class_hash.keys():
        print(ck,class_hash[ck], end=' ')
    print(' ')

    # how many of each subclass?
    subclass_hash = frames.subclass.value_counts()
    for sk in subclass_hash.keys():
        print(sk,subclass_hash[sk], end=' ')
    print(' ')

    # how many by each analyst?
    analyst_hash = frames.analyst.value_counts()
    for ak in analyst_hash.keys():
        print(ak,analyst_hash[ak], end=' ')
    print(' ')

    print('***********************\n')

    # Now create the day summary dataframes
    # What we really want to do here is build a CSV file that contains date, and number of each mainclass and subclass for that date
    # use mainclass and subclass from above
    
    # create a new dataframe to do daily counts of eventtype
    cols = ['date']
    for ck in class_hash.keys():
        cols.append(ck)
    for sk in subclass_hash.keys():
        cols.append(sk)
    #eventtype_df = pd.DataFrame(columns=['date','D','R','LV','r','e','l','h','t'])
    eventtype_df = pd.DataFrame(columns=cols)
    eventtype_df.set_index('date',inplace=True)

    # create a new dataframe to do daily counts by analyst
    cols2 = ['date']
    for ak in analyst_hash.keys():
        cols2.append(ak)
    analyst_df = pd.DataFrame(columns=cols2)
    analyst_df.set_index('date',inplace=True)

    for csvfile in list_of_csv_files:

        # read CSV into dataframe
        df = pd.read_csv(csvfile)
        # simplify column names to lowercase and no space
        df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')

        if len(df.mainclass.unique()): # number of mainclasses must be >0

            yyyy = int(csvfile[11:15])
            mm = int(csvfile[15:17])
            start_date = dt.date(yyyy, mm, 1)
            end_date = dt.date(yyyy, mm+1, 1)
            delta = dt.timedelta(days=1)
            while start_date < end_date: # loop over all days in this month
                yyyymmdd = start_date.strftime("%Y%m%d")
                print(yyyymmdd)
                sfiledatetime = pd.to_datetime(df.datetime)
                daycat = df.loc[sfiledatetime.dt.date == start_date] # get rows just for this day
                dayclass_hash = daycat.mainclass.value_counts() # mainclasses for this day
                daysubclass_hash = daycat.subclass.value_counts() # subclasses for this day
                dayanalyst_hash = daycat.analyst.value_counts() # analysts for this day
                # now go through the mainclass and subclass hashes from the overall summary, and make into a new dataframe here
                for ck in class_hash.keys():
                    if ck in dayclass_hash:
                        eventtype_df.at[yyyymmdd, ck] = dayclass_hash[ck] 
                    else:
                        eventtype_df.at[yyyymmdd, ck] = 0
                for sk in subclass_hash.keys():
                    if sk in daysubclass_hash:
                        eventtype_df.at[yyyymmdd, sk] = daysubclass_hash[sk] 
                    else:
                        eventtype_df.at[yyyymmdd, sk] = 0
                for ak in analyst_hash.keys():
                    if ak in dayanalyst_hash:
                        analyst_df.at[yyyymmdd, ak] = dayanalyst_hash[ak] 
                    else:
                        analyst_df.at[yyyymmdd, ak] = 0
                
                start_date += delta # add a day
        #nrows, ncolumns = df.shape
    print(eventtype_df)
    print(analyst_df)
    return eventtype_df, analyst_df

def generate_monthly_wav_csv(mm):
# function generate_monthly_wav_csv(monthdirs)
    print("Reading %s" % mm)
    wavfileslist = sorted(glob.glob(os.path.join(mm, "*")))
    outfile = mm[-13:-8] + 'wavfiles' + mm[-7:-3] + mm[-2:] + '.csv'
    print("...Creating %s" % outfile)
    if os.path.exists(outfile):
        print("CSV file %s already exists. Not overwriting. If you want to overwrite then please delete from the command line" % outfile)
        return
    fptr = open(outfile,'w')
    fptr.write('datetime,duration,wavfilepath,sampling_rate,npts,traceNum,traceID,fixedID,quality_factor,snr,highval,lowval\n') 
    for thiswavfile in wavfileslist:
        if thiswavfile.find('.png')>-1:
            continue
        try:
            st = obspy.read(thiswavfile)
        except:
            print("Cannot read wavfile %s" % thiswavfile )
            continue
        tmp = thiswavfile.split('WAV')
        shortwavfile = 'WAV' + tmp[1] 
        basename = os.path.basename(shortwavfile)
        pos = basename.find('S.')
        #if pos==15: # e.g. 9610-25-1005-41
        #    yyyy = int(basename[0:2]) + 1900
        #elif pos==18:
        #    yyyy = int(basename[0:4])
        #mm=int(basename[pos-13:pos-11])
        #dd=int(basename[pos-10:pos-8])
        #hh=int(basename[pos-7:pos-5])
        #mi=int(basename[pos-5:pos-3])
        #ss=int(basename[pos-2:pos])
        #filetime = dt.datetime(yyyy,mm,dd,hh,mi,ss)
        tracenum = 0
        for tr in st:
            tracetime = tr.stats.starttime
            tr2, quality_factor, snr = qc.compute_metrics(tr)
            duration = tr.stats.npts / tr.stats.sampling_rate
            fptr.write("%s,%s,%8.3f,%3d,%6d,%02d,%s,%s,%d,%.2f,%.2f,%.2f\n" \
                %  (tracetime.datetime, \
                    shortwavfile, \
                    duration, \
                    tr.stats.sampling_rate, tr.stats.npts, \
                    tracenum, tr.id, tr2.id, \
                    quality_factor, snr[0], snr[1], snr[2] ))
            tracenum += 1
            
    fptr.close()


if __name__ == "__main__":
    pass
    # could read in Sfiles then 
    # export to ObsPy event object
    #e = s1.to_css()
    #print e

