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
r.read_sfile() # print lines from the file to screen
vars(r) # print attributes
print(r) # fails because there needs to be an origin. Default behaviour needed.
r.parse_sfile() # populates attributes
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
  init__(self,reafile): # constructor
  read_sfile(self):
  print_sfile(self,lines):
  write_sfile(self,lines):
  parse_sfile(self): # populate attributes from reafile
  set_origintime(self,ot):
  set_location(self,lat,lon,depth):
  set_model(self,model):
  add_external_mag(self,mag,type,agency,info):
  add_comment(self,comment):
  set_line1(self,ot,model,distance,latitude,longitude,depth,nsta,rms,mag,fixed=True):
  add_waveform(self,filename):
  preferred_magnitude(self,agency=None,code=None):
  maximum_magnitude(self):
  __str__(self):
  to_obspy(self): # read WAV file into ObsPy stream

And functions:
  lntype(line):
  mag_commentln(mag,type,agency,info):
  event_list(startdate,enddate): # return a Python list of all events in REA tree


"""
# Import PYTHON modules we need
import os
#import sys
from datetime import *
#import time
#import re


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

def event_list(startdate,enddate):
  event_list=[]
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
  

  def read_sfile(self):
    try:
      rea_fh=open(self.file,"r")
      raw=rea_fh.read()
      lines=raw.split('\n')
      rea_fh.close()
      return lines
    except IOError:
      print(("Error: The specified file does not exist - %s" % (self.file)))
      raise e

  def write_sfile(self,lines):
    try:
      rea_fh=open(self.file,'w')
      for line in lines:
        rea_fh.write("%s\n" % (line))
      rea_fh.close()
    except IOError:
      print(("Error: The specified file does not exist - %s" % (self.file)))
      raise e
  
  def print_sfile(self,lines):
    for line in lines:
      print(line)
    return True

  def parse_sfile(self):
    lines=self.read_sfile()
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
  
  def set_origintime(self,ot):
    #ot must be a datetime object
    lines=self.read_sfile()
    if lntype(lines[0]) == '1':
      dstr=" %4d %2d%2dF%2d%2d%5.1f" % (ot.year,ot.month,ot.day,ot.hour,ot.minute,ot.second+ot.microsecond*1e-6)
      sufix=lines[0][20:len(lines[0])]
      lines[0]="%s%s" % (dstr,sufix)
      if len(lines[0])==80:
        self.write_sfile(lines)
        self.debug_msg="Modified Type 1 Line: \n%s" % lines[0]
        return True
      else:
        self.error_msg="Line 1 has incorrect length"      
        return False
    else:
      self.error_msg="Something funny happened and line one is not the first line"
      return False

  def set_location(self,lat,lon,depth):
    lines=self.read_sfile()
    if lntype(lines[0]) == '1':
      dstr="%7.3f%8.3f%5.1fFF" % (lat,lon,depth)
      prefix=lines[0][0:23]
      sufix=lines[0][45:len(lines[0])]
      lines[0]="%s%s%s" % (prefix,dstr,sufix)
      if len(lines[0])==80:
        self.write_sfile(lines)
        self.debug_msg="Modified Type 1 Line: \n%s" % lines[0]
        return True
      else:
        self.error_msg="Line 1 has incorrect length"      
        return False
    else:
      self.error_msg="Something funny happened and line one is not the first line"
      return False  

  def set_model(self,model):
    lines=self.read_sfile()
    if lntype(lines[0]) == '1':
       prefix=lines[0][0:20]
       sufix=lines[0][21:len(lines[0])]
       lines[0]="%s%s%s" % (prefix,model,sufix)
       self.write_sfile(lines)
       self.debug_msg="Modified Type 1 Line: \n%s" % lines[0]
       return True
    else:
      self.error_msg="Something funny happened and line one is not the first line"
      return False

  def add_external_mag(self,mag,type,agency,info):
  # Support for using comment lines to contain external magnitudes
    lines=self.read_sfile()
    ln3=mag_commentln(mag,type,agency,info)
    for index, line in enumerate(lines):
      if lntype(line) == 'I':
        lines.insert(index+1,ln3)
        self.debug_msg="Added ExtMag Comment Line: \n%s" % (ln3)
    self.write_sfile(lines)

  def add_comment(self,comment):
  # Support for adding a comment to an sfile
  # comment must be 78 characters or less
    lines=self.read_sfile()
    ln3=' '+comment+' '*(78-len(comment))+'3'    
    for index, line in enumerate(lines):
      if lntype(line) == 'I':
        lines.insert(index+1,ln3)
        self.debug_msg="Added Comment Line: \n%s" % (ln3)
    self.write_sfile(lines)

  def set_line1(self,ot,model,distance,latitude,longitude,depth,nsta,rms,mag,fixed=True):
    """ This function should be used sparingly used to fix corrupted line 1"""
    lines=self.read_sfile()
    ln1=" %4d %02d%02dF%02d%02d% 4.1f%s%s %7.3f%8.3f%5.1fFFOGS%3d%4.1f%4.1f OGS                1" % (ot.year,ot.month,ot.day,ot.hour,ot.minute,(ot.second+ot.microsecond*1e-6),model,distance,latitude,longitude,depth,nsta,rms,mag)
    lines[0]=ln1
    self.write_sfile(lines)
    self.debug_msg="Modified Type 1 Line: \n%s" % lines[0]    
    return True    

  def add_waveform(self,filename):
    lines=self.read_sfile()
    ln6=" "+filename
    sp=(79-len(ln6))*" "
    ln6=ln6+sp+"6"
    for index, line in enumerate(lines):
      if lntype(line) == 'I':
        lines.insert(index+1,ln6)
        self.debug_msg="Added waveform file line:\n%s" % (ln6)
    self.write_sfile(lines)
    return True

  def preferred_magnitude(self,agency=None,code=None):
    # code is the magnitude type
    select_agency=agency
#   print select_agency
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

  def maximum_magnitude(self):
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

  def to_obspy(self):
    """
    returns an obspy  stream object
    """
    import obspy.core
    
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

# End of REA Class 
