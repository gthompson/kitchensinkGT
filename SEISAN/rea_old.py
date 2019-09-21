""" rea.py
  This contains classes and methods for accessing and manipulating
 data in a SEISAN catalog single REA file.

 Austin Holland
 University of Oklahoma
 2010
"""

"""
Modifications by GLenn THompson 2014
1. There were a mixture of tabs and spaces which meant my version of Python couldn't load the module.
Stripped out the tabs and replaced with spaces, then realigned code throughout whole file to prevent
identation errors.
2. 

Methods used are:
from rea import REA
r = REA('/raid/data/seisan/REA/MVOE_/2002/01/31-2356-00L.S200201') # create the object
r.read_file() # print lines from the file to screen
vars(r) # print attributes
print(r) # fails because there needs to be an origin. Default behaviour needed.
r.parse_reaheader() # populates attributes
r.to_obspy()
r.plot()
Got error:
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/raid/apps/antelope/python2.7.2/lib/python2.7/site-packages/obspy-0.8.4_705_gd3b2fe_dirty-py2.7.egg/obspy/core/stream.py", line 1083, in plot
    waveform = WaveformPlotting(stream=self, *args, **kwargs)
  File "/opt/antelope/python2.7.2/lib/python2.7/site-packages/obspy.imaging-0.7.0-py2.7.egg/obspy/imaging/waveform.py", line 59, in __init__
    raise IndexError(msg)
IndexError: Empty stream object
>>> print st
0 Trace(s) in Stream:

r.to_obspy() obviously is not working!

Here is a full list of REA class methods:
  __init__(self,reafile): # constructor
  read_file(self):
  print_file(self,lines):
  write_file(self,lines):
  parse_reaheader(self): # populate attributes from reafile
  fix_depth(self,depth):
  free_depth(self):
  fix_origintime(self,ot):
  fix_location(self,lat,lon,depth):
  set_velmod(self,model):
  external_mag(self,mag,type,agency,info):
  add_comment(self,comment):
  gen_ln1(self,ot,model,distance,latitude,longitude,depth,nsta,rms,mag,fixed=True):
  add_waveform(self,filename):
  phase_traveltime(self,h,m,s):
  parse_Amb(self):
  parse_AML(self):
  parse_pmotions(self):
  parse_phasetimes(self):
  calculate_mblg(self,quiet=False):
  mblg_recalc(self):
  preferred_magnitude(self,agency=None,code=None,force_mblg_recalc=False):
  maximum_magnitude(self,force_mblg_recalc=False):
  __str__(self):
  catalog_report(self):
  to_simpleList(self):
  to_GEORSS(self):
  to_EQXML(self,action=None):
  to_CSVsimple(self):
  to_XMLCatalog(self):
  to_HTML(self,filename=None):
  to_obspy(self): # read WAV file into ObsPy stream
  parse_traveltimes(self,phase):
  to_SimpleEQ(self):

And functions:
  lntype(line):
  mag_commentln(mag,type,agency,info):
  sec2isecmsec(seconds):
  velocity(dist,time):
  mean(array):
  seidb_traverse(readir,dbname,type):
  event_list(startdate,enddate): # return a Python list of all events in REA tree
  eqxml_isoformat(dt):
  eqxml_magkey(type):
  get_eqxmlID():


"""
# Import PYTHON modules we need
import os
import sys
from datetime import *
import time
import eqcatalog
import math
import re
#from eqdb.py import *


#Useful global variables
magtype_map={"L":"ML","b":"mb","B":"mB","s":"Ms","S":"MS","W":"MW","G":"mbLg","C":"Md"}

def lntype(line):
  if len(line) == 80:
    type=line[79]
  else:
    type=''
  return type

def mag_commentln(mag,type,agency,info):
  str=" ExtMag %4.1f%1s%3s%63s3" % (mag,type,agency,info)
  return str

def sec2isecmsec(seconds):
  isec=int(seconds)
  msec=int((seconds-isec)*1e6)
  return isec,msec

def velocity(dist,time):
  return dist/time       

def mean(array):
  n=len(array)
  if n > 0:
    return (math.fsum(array)/n)
  else:
    return 0.0

def seidb_traverse(readir,dbname,type):
  """full_filenames=seidb_traverse(readir,dbname,type)
  This will return a list of fully qualified path and filenames for the SEISAN database
  identified by the dbname and type is the SEISAN location code L,R, or D."""
  _sdir=""
  _sstr="%s.S" % (type)
  flist=[]
  if readir[-1]=='/':
    _sdir=readir+dbname
  else:
    _sdir=readir+'/'+dbname
  for root, dirs, files in os.walk(_sdir):
    for filename in files:
      if re.search(_sstr,filename):
        flist.append(root+'/'+filename)
  return flist

def event_list(startdate,enddate):
  event_list=[]
  #reapath="/home/analyst/REA/OGS__/"
  reapath="/raid/data/seisan/REA/MVOE_/"
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
      dir="%s%d/%02d/" % (reapath,year,month)
      #print dir
      flist=glob(dir+"*L.S*")
      event_list.extend(flist)
  return event_list  

class REA:

  def __init__(self,reafile):
    if os.path.exists(reafile):
      self.file=reafile
      self.file_mtime=os.path.getmtime(reafile) 
      self.wave_files=[]
      self.error_msg=""
      self.evtype={}
      self.magnitude=[]
      self.magnitude_type=[]
      self.magnitude_agency=[]
      self.line1_proc=False
      self.error={}
      self.url=''
      self.ha={}
      self.focmec={}
      self.maximum_intensity=0
      self.debug_msg="Initialized instance of REA Class"
      self.error_msg=''
      self.origin_time = None
    else:
      self.error_msg="REA file does not exist"
  

  def read_file(self):
    try:
      rea_fh=open(self.file,"r")
      raw=rea_fh.read()
      lines=raw.split('\n')
      rea_fh.close()
      return lines
    except IOError:
      print("Error: The specified file does not exist - %s" % (self.file))
      raise e

  def print_file(self,lines):
    for line in lines:
      print(line)
    return True

  def write_file(self,lines):
    try:
      rea_fh=open(self.file,'w')
      for line in lines:
        rea_fh.write("%s\n" % (line))
      rea_fh.close()
    except IOError:
      print("Error: The specified file does not exist - %s" % (self.file))
      raise e
  
  def parse_reaheader(self):
    lines=self.read_file()
    for line in lines:
      if lntype(line) != '':
        # Process Type 1 line
        if lntype(line) == '1':
          if not self.line1_proc:
            self.line1_proc=True
            year=int(line[1:5])
            month=int(line[6:8])
            day=int(line[8:10])
            hour=int(line[11:13])
            min=int(line[13:15])
            isec=int(line[16:18])
            if isec == 60: # 
               min+=1
               isec=isec-60
            self.origin_time=datetime(year,month,day,hour,min,isec)
            self.evtype['dist']=line[21]
            self.velmod=line[20]
            self.evtype['qual']=line[22]
            if line[23:30]!="       ":
              self.latitude=float(line[23:30])
            if line[30:38]!="        ":
              self.longitude=float(line[30:38])
            if line[38:43]!="     ":
              self.depth=float(line[38:43])
            self.z_indicator=line[43]
            self.agency=line[45:48]
            if line[48:51]!="   ":
              self.no_sta=int(line[48:51])
            else:
              self.no_sta=0
            if line[51:55]!="    ":
              self.rms=float(line[51:55])
            if line[55:59] != '    ':
              self.magnitude.append(float(line[55:59]))
              self.magnitude_type.append(magtype_map[line[59]])
              self.magnitude_agency.append(line[60:63])
            if line[63:67] != '    ':
              self.magnitude.append(float(line[63:67]))
              self.magnitude_type.append(magtype_map[line[67]])
              self.magnitude_agency.append(line[68:71])
            if line[71:75] != '    ':
              self.magnitude.append(float(line[71:75]))
              self.magnitude_type.append(magtype_map[line[75]])
              self.magnitude_agency.append(line[76:79])
          else: #If there is an additional line1 process the additional magnitudes
            if line[55:59] != '    ':
              self.magnitude.append(float(line[55:59]))
              self.magnitude_type.append(magtype_map[line[59]])
              self.magnitude_agency.append(line[60:63])
            if line[63:67] != '    ':
              self.magnitude.append(float(line[63:67]))
              self.magnitude_type.append(magtype_map[line[67]])
              self.magnitude_agency.append(line[68:71])
            if line[71:75] != '    ':
              self.magnitude.append(float(line[71:75]))
              self.magnitude_type.append(magtype_map[line[75]])
              self.magnitude_agency.append(line[76:79])

        # Process Type I line, ID Line
        if lntype(line) == 'I':
          self.last_action=line[8:11]
          self.action_time=line[12:26]
          self.analyst=line[30:32]
          self.id=int(line[60:74])

        # Process Type 3 line, Comment Lines for extra information
        if lntype(line) == '3':
          if line[1:7]== 'ExtMag':
            self.magnitude.append(float(line[8:12]))
            self.magnitude_type.append(magtype_map[line[12]])
            self.magnitude_agency.append(line[13:16])
          if line[1:4]=='URL':
            self.url=line[5:78].rstrip()

        # Process Type E line, Hyp error estimates
        if lntype(line) == 'E':
          #print line
          self.gap=int(line[5:8])
          self.error['origintime']=float(line[14:20])
          self.error['latitude']=float(line[24:30].strip())
          self.error['longitude']=float(line[32:38].strip())
          self.error['depth']=float(line[38:43].strip())
          self.error['covxy']=float(line[43:55].strip())
          self.error['covxz']=float(line[55:67].strip())
          self.error['covyz']=float(line[67:79].strip())

        # Process Type H line, High accuracy line
        if lntype(line) == 'H':
          self.ha['ot_sec']=float(line[16:22])
          self.origin_time=datetime(int(line[1:5]),int(line[6:8]),int(line[8:10]),
            int(line[11:13]),int(line[13:15]),int(self.ha['ot_sec']),
            int((self.ha['ot_sec']-int(self.ha['ot_sec']))*1.e6))
          self.ha['latitude']=float(line[23:32].strip())
          self.ha['longitude']=float(line[33:43].strip())
          self.ha['depth']=float(line[44:52].strip())
          self.ha['rms']=float(line[53:59].strip())
          self.latitude=self.ha['latitude']
          self.longitude=self.ha['longitude']
          self.depth=self.ha['depth']
          self.rms=self.ha['rms']

        # Process Type F line, Fault plane solution
        # Format has changed need to fix AAH - 2011-06-23
        if lntype(line) == 'F' and 'dip' not in self.focmec:
          self.focmec['strike']=float(line[0:10])
          self.focmec['dip']=float(line[10:20])
          self.focmec['rake']=float(line[20:30])
          #self.focmec['bad_pols']=int(line[60:66])
          self.focmec['agency']=line[66:69]
          self.focmec['source']=line[70:77]
          self.focmec['quality']=line[77]
  
        # Process Type 6 lines, Waveform file lists
        # Build an array of waveform files
        if lntype(line) == '6':
          lvals=line[1:79].split()
          for val in lvals:
            self.wave_files.append(val)

        # Process Type 2 line, Macroseismic Intensity Information
        if lntype(line) == '2':
          self.maximum_intensity=int(line[27:29])
 
#        line=rea_fh.readline()
#        print line[79]

    return True
  
  def fix_depth(self,depth):
    lines=self.read_file()
    if lntype(lines[0]) == '1':
       dstr="%5.1fF" % (depth)
       prefix=lines[0][0:38]
       sufix=lines[0][44:len(lines[0])]
       lines[0]="%s%s%s" % (prefix,dstr,sufix)
       self.write_file(lines)
       self.debug_msg="Modified Type 1 Line: \n%s" % lines[0]
       return True
    else:
      self.error_msg="Something funny happened and line one is not the first line"
      return False
      
  def free_depth(self):
    if self.z_indicator=='F':
      lines=self.read_file()
      if lntype(lines[0]) == '1':
        dstr="%5.1f " % (self.depth)
        prefix=lines[0][0:38]
        sufix=lines[0][44:len(lines[0])]
        lines[0]="%s%s%s" % (prefix,dstr,sufix)
        self.write_file(lines)
        self.debug_msg="Modified Type 1 Line: \n%s" % lines[0]
        return True
      else:
        self.error_msg="Something funny happened and line one is not the first line"
        return False   
    else:
      return False   

  def fix_origintime(self,ot):
    #ot must be a datetime object
    lines=self.read_file()
    if lntype(lines[0]) == '1':
      dstr=" %4d %2d%2dF%2d%2d%5.1f" % (ot.year,ot.month,ot.day,ot.hour,ot.minute,ot.second+ot.microsecond*1e-6)
      sufix=lines[0][20:len(lines[0])]
      lines[0]="%s%s" % (dstr,sufix)
      if len(lines[0])==80:
        self.write_file(lines)
        self.debug_msg="Modified Type 1 Line: \n%s" % lines[0]
        return True
      else:
        self.error_msg="Line 1 has incorrect length"      
        return False
    else:
      self.error_msg="Something funny happened and line one is not the first line"
      return False

  def fix_location(self,lat,lon,depth):
    lines=self.read_file()
    if lntype(lines[0]) == '1':
      dstr="%7.3f%8.3f%5.1fFF" % (lat,lon,depth)
      prefix=lines[0][0:23]
      sufix=lines[0][45:len(lines[0])]
      lines[0]="%s%s%s" % (prefix,dstr,sufix)
      if len(lines[0])==80:
        self.write_file(lines)
        self.debug_msg="Modified Type 1 Line: \n%s" % lines[0]
        return True
      else:
        self.error_msg="Line 1 has incorrect length"      
        return False
    else:
      self.error_msg="Something funny happened and line one is not the first line"
      return False  

  def set_velmod(self,model):
    lines=self.read_file()
    if lntype(lines[0]) == '1':
       prefix=lines[0][0:20]
       sufix=lines[0][21:len(lines[0])]
       lines[0]="%s%s%s" % (prefix,model,sufix)
       self.write_file(lines)
       self.debug_msg="Modified Type 1 Line: \n%s" % lines[0]
       return True
    else:
      self.error_msg="Something funny happened and line one is not the first line"
      return False

  def external_mag(self,mag,type,agency,info):
  # Support for using comment lines to contain external magnitudes
    lines=self.read_file()
    ln3=mag_commentln(mag,type,agency,info)
    for index, line in enumerate(lines):
      if lntype(line) == 'I':
        lines.insert(index+1,ln3)
        self.debug_msg="Added ExtMag Comment Line: \n%s" % (ln3)
    self.write_file(lines)

  def add_comment(self,comment):
  # Support for adding a comment to an sfile
  # comment must be 78 characters or less
    lines=self.read_file()
    ln3=' '+comment+' '*(78-len(comment))+'3'    
    for index, line in enumerate(lines):
      if lntype(line) == 'I':
        lines.insert(index+1,ln3)
        self.debug_msg="Added Comment Line: \n%s" % (ln3)
    self.write_file(lines)

  def gen_ln1(self,ot,model,distance,latitude,longitude,depth,nsta,rms,mag,fixed=True):
    """ This function should be used sparingly used to fix corrupted line 1"""
    lines=self.read_file()
    ln1=" %4d %02d%02dF%02d%02d% 4.1f%s%s %7.3f%8.3f%5.1fFFOGS%3d%4.1f%4.1f OGS                1" % (ot.year,ot.month,ot.day,ot.hour,ot.minute,(ot.second+ot.microsecond*1e-6),model,distance,latitude,longitude,depth,nsta,rms,mag)
    lines[0]=ln1
    self.write_file(lines)
    self.debug_msg="Modified Type 1 Line: \n%s" % lines[0]    
    return True    


  def add_waveform(self,filename):
    lines=self.read_file()
    ln6=" "+filename
    sp=(79-len(ln6))*" "
    ln6=ln6+sp+"6"
    for index, line in enumerate(lines):
      if lntype(line) == 'I':
        lines.insert(index+1,ln6)
        self.debug_msg="Added waveform file line:\n%s" % (ln6)
    self.write_file(lines)
    return True
    
  def phase_traveltime(self,h,m,s):
    day=0
    if h >= 24:
      day += 1
      h -= 24
  # This may have to be modified some, but appears to work
    osec,omsec=sec2isecmsec(self.ha['ot_sec'])
    haot=datetime(self.origin_time.year,self.origin_time.month,self.origin_time.day,self.origin_time.hour,self.origin_time.minute,osec,omsec)
    psec,pmsec=sec2isecmsec(s)
    pt=datetime(self.origin_time.year,self.origin_time.month,self.origin_time.day+day,h,m,psec,pmsec)
    tt=pt - haot
    ttseconds=tt.seconds + (tt.microseconds/1e6)
    return ttseconds
  

  def parse_Amb(self):
    lines=self.read_file()
    amps=[]
    for line in lines:
      if lntype(line) == ' ' or lntype(line) == '4':
        if line[10:14] == 'IAmb':
          tmp={}
          #print line[1:6]
          tmp['sta']=line[1:6]
          tmp['inst']=line[6]
          tmp['comp']=line[7]
          tmp['hour']=int(line[18:20])
          tmp['min']=int(line[20:22])
          tmp['sec']=float(line[22:28])
#          print line[33:40]
          try:
            tmp['amplitude']=float(line[33:40])
            tmp['period']=float(line[41:45])
            tmp['distance']=float(line[70:75])
            tmp['azimuth']=int(line[76:79])
            amps.append(tmp)
          except:
            print("ERROR READING LINE")
            print(line)
            tmp={}
    return amps

  def parse_AML(self):
    lines=self.read_file()
    amps=[]
    for line in lines:
      if lntype(line) == ' ' or lntype(line) == '4':
        if line[10:14] == 'IAML':
          if line[33:40]!='       ':
            tmp={}
            #print line[1:6]
            tmp['sta']=line[1:6]
            tmp['inst']=line[6]
            tmp['comp']=line[7]
            tmp['hour']=int(line[18:20])
            tmp['min']=int(line[20:22])
            tmp['sec']=float(line[22:28])
#           print line[33:40]
            tmp['amplitude']=float(line[33:40])
            tmp['period']=float(line[41:45])
            tmp['distance']=float(line[70:75])
            tmp['azimuth']=int(line[76:79])
            amps.append(tmp)
    return amps


  def parse_pmotions(self):
    lines=self.read_file()
    mot=[]
    for line in lines:
      if lntype(line) == ' ' or lntype(line) == '4':
        if line[10] == 'P':
          if line[16]!=' ':
            tmp={}
            #print line[1:6]
            tmp['sta']=line[1:6]
            tmp['qual']=line[9]
            tmp['mot']=line[16]
            tmp['azm']=float(line[76:79])
            tmp['ain']=float(line[56:60])
            mot.append(tmp)
    return mot

  def parse_phasetimes(self):
    origintime=self.origin_time
    lines=self.read_file()
    tmp={}
    for line in lines:
      if lntype(line) == ' ' or lntype(line) == '4':
        sta=line[1:6].strip()
        
        if line[10] == 'P':
          if not sta in list(tmp.keys()) and len(sta)>2:
                tmp[sta]={}
          try:
            tmp[sta]['dist']=float(line[71:75].strip())
            tmp[sta]['pqual']=line[9]
            tmp[sta]['pmot']=line[16]
            tmp[sta]['pdesc']=line[11]
            hr,mn,sec=int(line[18:20]),int(line[20:22]),float(line[22:28])
            if hr==24:
              hr=0
              day=origintime.day+1
            else:
              day=origintime.day
            tmp[sta]['ptime']=datetime(origintime.year,origintime.month,
                day,hr,mn,int(sec),int((sec-int(sec))*1e6))
          except:
            pass
        if line[10] == 'S':
          if not sta in list(tmp.keys()) and len(sta)>2:
                tmp[sta]={}
          try:
            tmp[sta]['dist']=float(line[71:75].strip())
            tmp[sta]['squal']=line[9]
            tmp[sta]['smot']=line[16]
            tmp[sta]['sdesc']=line[11]
            hr,mn,sec=int(line[18:20]),int(line[20:22]),float(line[22:28])
            if hr==24:
              hr=0
              day=origintime.day+1
            else:
              day=origintime.day
            tmp[sta]['stime']=datetime(origintime.year,origintime.month,
                day,hr,mn,int(sec),int((sec-int(sec))*1e6))
          except:
            pass

    return tmp

        
  def calculate_mblg(self,quiet=False):
    Q=1000.0 # Default Q value
    amps=self.parse_Amb()
    mags=[]
    if not quiet:
      print("mbLg Magnitude Calculation summary")
      print("\t Sta \t Dist\t Vel \tAmp(um)\t Per\t mbLg\tComments")
    for amp in amps:
      comment=""
      lgtt=self.phase_traveltime(amp['hour'],amp['min'],amp['sec'])
      lgvel=velocity(amp['distance'],lgtt)
      gamma=math.pi/(Q*lgvel*amp['period'])
      amplitude=(amp['amplitude']/2)/1e3 # Convert from nanometers to micrometers
      mblg=2.96+0.833*math.log10(amp['distance']/10.0)+.4343*gamma*amp['distance']+math.log10(amplitude)
      if amp['distance'] >= 55.6 and amp['distance'] <= 1000.0:
        if lgvel >= 3.0 and lgvel <= 3.7:
          if amp['period'] >= 0.33 and amp['period'] <= 1.333: # Set our cut off to a full octave above and below 1.5 Hz of the WWSN SP peak response
             mags.append(mblg)
          else:
             comment+="Bad Period, "
        else:
          comment+="Bad Velocity, "
      else:
        comment+="Bad Distance, "
      if not quiet:
        print("\t%s\t%5.1f\t%5.2f\t%7.5f\t%4.2f\t%5.2f\t%s" %(amp['sta'],amp['distance'],lgvel,amplitude,amp['period'],mblg,comment))
    if not quiet:
      print("\nmbLg= %4.1f" % (mean(mags)))
    return mags

  def mblg_recalc(self):
  #This function was added to make sure that the mbLg reported is determined from the most
  #current data and not what is necessarily stored in the SEISAN
    mblg_index=-1
    mags=[]
    for i,magnitude in enumerate(self.magnitude):
      if self.magnitude_type[i]=="mbLg" and self.magnitude_agency[i]=="OGS":
        mblg_index=i
    if mblg_index>=0:
      mags=self.calculate_mblg(quiet=True)
      self.magnitude[mblg_index]=mean(mags)
    else:
      if len(mags)>0:
        self.magnitude.append(mean(mags))
        self.magnitude_type.append("mbLg")
        self.magnitude_agency("OGS")
    return True

  def preferred_magnitude(self,agency=None,code=None,force_mblg_recalc=False):
    # code is the magnitude type
    select_agency=agency
#   print select_agency
    if force_mblg_recalc:
      self.mblg_recalc()
    agency_equiv={' GS':' US'}
    pref_order=[
        {'agency':' US','type':'MW'},
        {'agency':' US','type':'MS'},
        {'agency':'OGS','type':'ML'},
        {'agency':' US','type':'mbLg'},
        {'agency':'OGS','type':'mbLg'},
        {'agency':'OGS','type':'Md'},
        {'agency':' US','type':'Md'}]
    mag=0.0
    mtype="NA"
    agency="NA"
    for j,magid in enumerate(pref_order):
      for i,magnitude in enumerate(self.magnitude):
        for key in list(agency_equiv.keys()):
          if self.magnitude_agency[i]==key:
            self.magnitude_agency[i]=agency_equiv[key]
        if magid['agency']==self.magnitude_agency[i] and magid['type']==self.magnitude_type[i]:
          if select_agency:
            if (agency=="NA") and (self.magnitude_agency[i]==select_agency):
              #Only update values if we haven't yet so the order of presedence is preserved
              if code:
                if code==self.magnitude_type[i]:
                  mag=magnitude
                  mtype=self.magnitude_type[i]
                  agency=self.magnitude_agency[i]
              else:
                mag=magnitude
                mtype=self.magnitude_type[i]
                agency=self.magnitude_agency[i]

          else:
            if agency=="NA":
              #Only update values if we haven't yet so the order of presedence is preserved
              if code:
                if code==self.magnitude_type[i]:
                  mag=magnitude
                  mtype=self.magnitude_type[i]
                  agency=self.magnitude_agency[i]
                else:
                  mag=magnitude
                  mtype=self.magnitude_type[i]
                  agency=self.magnitude_agency[i]
    return mag,mtype,agency


  def maximum_magnitude(self,force_mblg_recalc=False):
    if force_mblg_recalc:
      self.mblg_recalc()
    mag=0.0
    mtype="NA"
    agency="NA"
    for i,magnitude in enumerate(self.magnitude):
      if magnitude > mag:
        mag=magnitude
        mtype=self.magnitude_type[i]
        agency=self.magnitude_agency[i]

    return mag,mtype,agency    
      

  def __str__(self):
    # GTHO 20140319 Added test here to see if origin_time defined
    if self.origin_time:
      return "%s (UTC) - Magnitude %.1f" % (self.origin_time.ctime(),self.preferred_magnitude()[0])
    else:
      return "No origin time"

  def to_EQXML(self,action=None):
    self.parse_reaheader()
    """This function converts a seisan REA file to an eqxml document."""
    from lxml import etree
    
    root=etree.Element("EQMessage",xmlns="http://www.usgs.gov/ansseqmsg")
    srccode=etree.SubElement(root,"Source")
    srccode.text="OK"
    sender=etree.SubElement(root,"Sender")
    sender.text="Oklahoma Geological Survey (OGS)"
    msgsrc=etree.SubElement(root,"MsgSrc")
    msgsrc.text="OGS"
    msgid=etree.SubElement(root,"MsgIdent")
    eqxml_id=get_eqxmlID()
    msgid.text=str(eqxml_id)
    sent=etree.SubElement(root,"Sent")
    sent.text=eqxml_isoformat(datetime.now())
    event=etree.SubElement(root,"Event")
    evsrc=etree.SubElement(event,"DataSource")
    evsrc.text="OK"
    evid=etree.SubElement(event,"EventID")
    evid.text=str(self.id)
    #if action!=None:
    #  evaction=etree.SubElement(event,"Action")
    #  evaction.text=action
    type=etree.SubElement(event,"Type")
    type.text="Earthquake"
    scope=etree.SubElement(event,"Scope")
    scope.text="Internal"
    origin=etree.SubElement(event,"Origin")
    ot=etree.SubElement(origin,"Time")
    ot.text=eqxml_isoformat(self.origin_time)
    lat=etree.SubElement(origin,"Latitude")
    lat.text=str(self.latitude)
    lon=etree.SubElement(origin,"Longitude")
    lon.text=str(self.longitude)
    depth=etree.SubElement(origin,"Depth")
    depth.text=str(self.depth)
    rms=etree.SubElement(origin,"StdError")
    rms.text=str(self.rms)
    gap=etree.SubElement(origin,"AzimGap")
    gap.text=str(self.gap)
    sta=etree.SubElement(origin,"NumStaUsed")
    sta.text=str(self.no_sta)
    oterr=etree.SubElement(origin,"OTError")
    oterr.text=str(self.error['origintime'])
    laterr=etree.SubElement(origin,"LatError")
    laterr.text=str(self.error['latitude'])
    lonerr=etree.SubElement(origin,"LonError")
    lonerr.text=str(self.error['longitude'])
    if self.error['depth']>0.0:
      zerr=etree.SubElement(origin,"DepthError")
      zerr.text=str(self.error['depth'])
      method=etree.SubElement(origin,"DepthMethod")
      method.text="Free"
    else :
      method=etree.SubElement(origin,"DepthMethod")
      method.text="Fixed"
    status=etree.SubElement(origin,'Status')
    status.text="Reviewed"
    mag={}
    msrc={}
    mtype={}
    mval={}
    for i in range(0,len(self.magnitude)):
      mag[i]=etree.SubElement(origin,"Magnitude")
      msrc[i]=etree.SubElement(mag[i],"SourceKey")
      msrc[i].text="OK"
      mtype[i]=etree.SubElement(mag[i],"TypeKey")
      mtype[i].text=eqxml_magkey(self.magnitude_type[i])
      mval[i]=etree.SubElement(mag[i],"Value")
      mval[i].text=str(self.magnitude[i])
    #Add phase picks
    
    #Write to file
    fname="OGS_%04d_EQ.XML" % (eqxml_id)
    fh=open('/home/analyst/etc/eqxml/'+fname,'w')
    fh.write(etree.tostring(root, pretty_print=True))
    fh.close()
    return eqxml_id

  def to_CSVsimple(self):
    str="%s,%f,%f,%f,%f\n" % (self.origin_time.isoformat("T"),self.ha['longitude'],self.ha['latitude'],self.ha['depth'],self.preferred_magnitude()[0])
    return str
      
  def to_XMLCatalog(self):
    """ This function creates an xml eqcatalog entry for an rea file and returns the xmlstr value"""
    if self.catalog_report():
      mag,type,agency=self.maximum_magnitude()
      cat=dict(entry='EQ',id=self.id)
      cat['latitude']="%.4f" % (self.ha['latitude'])
      cat['longitude']="%.4f" % (self.ha['longitude'])
      cat['depth']="%.2f" % (self.ha['depth'])
      cat['magnitude']="%.1f" % (mag)
      cat['magnitudetype']=type
      # GTHO 20140318: add checks for existence of self.error keys as getting error here
      if 'latitude' in self.error:
        cat['laterror']="%.2f" % (self.error['latitude'])
      if 'longitude' in self.error:
        cat['lonerror']="%.2f" % (self.error['longitude'])
      if 'depth' in self.error:
        cat['deptherror']="%.2f" % (self.error['depth'])
      if 'origintime' in self.error:
        cat['oterror']="%.2f" % (self.error['origintime'])
      cat['origintime']="%s" % (self.origin_time.isoformat(" "))
      if eqcatalog.isdst(self.origin_time):
        ltime= self.origin_time - timedelta(hours=5.0)
        cat['localtime']="%s CDT" % (ltime.isoformat(" "))
      else:
        ltime= self.origin_time - timedelta(hours=6.0)
        cat['localtime']="%s CST" % (ltime.isoformat(" "))
      cat['nstas']="%i" % (self.no_sta)
      cat['rms']="%.2f" % (self.ha['rms'])
      cat['MaxMMI']=eqcatalog.roman_intensity(self.maximum_intensity)
      cat['updated']=datetime.utcfromtimestamp(self.file_mtime).isoformat("T")
      return eqcatalog.dict2xml(cat,2)
    else:
      return ''
    
  def to_HTML(self,filename=None):
    """ object.to_HTML(filename)
    If filename is not defined then the output is directed to standard out.
    """
    from eqcatalog import calc_localtime,isdst
    if not filename:
      fh=sys.stdout
    else:
      fh=open(filename,'w')
    fh.write("<html>\n")
    fh.write("<div style=\"width:600px,padding:20px;\">\n")
    fh.write("<h2>Earthquake Details: Magnitude %.1f - %s</h2>\n" % (self.preferred_magnitude()[0],self.origin_time.isoformat(" ")))
    fh.write("<table style=\"border: 1px solid black;\">\n")
  
    if isdst(self.origin_time):
      ltimestr="%s (CDT)" % (calc_localtime(self.origin_time))
    else:
      ltimestr="%s (CST)" % (calc_localtime(self.origin_time))

    fh.write("<tr><td style=\"font-weight:bold;vertical-align:top;\">Origin Time</td><td>%s (UTC)<br />%s<br />Uncertainty: %.3f seconds</td>\n" % (self.origin_time.isoformat(" "),ltimestr,self.error['origintime']))
    fh.write("<tr><td style=\"font-weight:bold;vertical-align:top;\">Location</td><td>Latitude: %.4f&deg; &plusmn; %.2f (km)<br />Longitude: %.4f&deg; &plusmn; %.2f (km) <br />Depth: %.2f &plusmn; %.2f km</td>\n" % (self.latitude,self.error['latitude'],self.longitude,self.error['longitude'],self.depth,self.error['depth']))
    magstr=''
    for i in range(0,len(self.magnitude)): # Build magnitude string
      magstr+="%.1f %s %s<br />" % (self.magnitude[i],self.magnitude_type[i],self.magnitude_agency[i])
    fh.write("<tr><td style=\"font-weight:bold;vertical-align:top;\">Magnitude</td><td>%s</td>\n" % (magstr))
    fh.write("<tr><td style=\"font-weight:bold;vertical-align:top;\">Maximum Modified Mercalli Intensity</td><td>%s</td>\n" % (eqcatalog.roman_intensity(self.maximum_intensity)))
    if 'strike' in self.focmec:
      mechstr="Strike: %.1f&deg;<br />Dip: %.1f&deg;<br />Rake: %.1f&deg;<br />Agency: %s<br />Method: %s,<br />" % (self.focmec['strike'],self.focmec['dip'],self.focmec['rake'],self.focmec['agency'],self.focmec['source'])
      fh.write("<tr><td style=\"font-weight:bold;vertical-align:top;\">Focal Mechanism</td><td>%s</td>\n" % (mechstr))
    fh.write("<tr><td style=\"font-weight:bold;vertical-align:top;\">Number of Stations</td><td>%d</td>\n" % (self.no_sta))
    fh.write("<tr><td style=\"font-weight:bold;vertical-align:top;\">Location RMS</td><td>%.2f</td>\n" % (self.rms))
    fh.write("<tr><td style=\"font-weight:bold;vertical-align:top;\">Last Updated</td><td>%s UTC</td>\n" % (datetime.utcfromtimestamp(self.file_mtime).isoformat(" ")))

    fh.write("</table>")

    #Write out the end of file and close
    fh.write("</div>\n</html>\n")
    if filename:
      fh.close()
    return True

  def to_obspy(self):
    """
    returns an obspy  stream object
    """
    import obspy.core
    
    #dir_top = '/home/analyst/WAV/'
    dir_top = '/raid/data/seisan/WAV/MVOE_/'
    years = os.walk(dir_top).next()[1]
    years.append('')
    st = obspy.core.Stream()
    for wave in self.wave_files:
      for year in years:
        fullyear = os.path.join(dir_top, year)
        #print fullyear
        months = os.walk(fullyear).next()[1]
        months.append('')
        for month in months:
          if len(month)>2:
            continue
          fullmonth = os.path.join(fullyear, month)
          #print fullmonth
          fullwave = os.path.join(fullmonth, wave)
          if os.path.exists(fullwave):
            st += obspy.read(fullwave)
            print("local pathname stored")
            continue
          elif os.path.exists(wave):
            st += obspy.read(wave)
            print("full pathname stored")
            continue
    return st

  def parse_traveltimes(self,phase):
    """
    Creates a dictionary of traveltimes
    of the format:
       {'receiver A name': travel time in seconds,
        'receiver B name': travel time in seconds}

    takes argument of self and phase (p,P,s,S)
    """
    
    self.parse_reaheader()
    phs=self.parse_phasetimes(self.origin_time)
    tts={}
    for sta in phs:
      if phase == 'p' or phase == 'P':
        try:
          sec=(phs[sta]['ptime']-self.origin_time).seconds
          msecs=float((phs[sta]['ptime']-self.origin_time).microseconds)/1000000
          tts[sta]=sec + msecs
        except:
          pass
      elif phase == 's' or phase == 'S':
        try:
          sec=(phs[sta]['stime']-self.origin_time).seconds
          msecs=float((phs[sta]['stime']-self.origin_time).microseconds)/1000000
          tts[sta]=sec + msecs
        except:
          pass
    return tts

"""
  def to_SimpleEQ(self):
    ret = SimpleEQ()
    ret.id=self.id
    ret.rea_file=self.file
    ret.last_mod=self.file_mtime
    ret.origin_time=self.origin_time
    ret.latitude=self.ha['latitude']
    ret.longitude=self.ha['longitude']
    ret.depth=self.ha['depth']
    ret.rms=self.ha['rms']
    ret.magnitude=self.maxmimum_magnitude()
    ret.err_lat=self.error['latitude']
    ret.err_lon=self.error['longitude']
    ret.err_depth=self.error['depth']
    ret.err_ot=self.error['origintime']
    ret.nstas=self.no_sta
    ret.intesity=self.maximum_intensity
    # Need to implement county determination
    ret.county="NA"

    return ret"""

# End of REA Class and define more local methods 
def eqxml_isoformat(dt):
  """ eqxml_tformat(datetime_instance)
  takes a datetime object and makes sure it gets represented properly for eqxml"""
  str="%d-%02d-%02dT%02d:%02d:%06.3fZ" % (dt.year,dt.month,dt.day,dt.hour,dt.minute,
    (dt.second+dt.microsecond/1.e6))
  return str
    
def eqxml_magkey(type):
  type_map={"ML":"Ml","mb":"mb","Ms":"Ms","MW":"Mw","mbLg":"mblg","Md":"Md"}
  return type_map[type]
  
def get_eqxmlID():
  fname='/home/analyst/etc/eqxml_evid.txt'
  idfile=open(fname,'r')
  id=idfile.readline()
  idfile.close()
  id=int(id.strip())
  id+=1
  
  f=open(fname,'w')
  f.write("%i\n" % (id))
  f.close()
  return id
