#! /usr/bin/env python
""" Converttime functions, Glenn Thompson, 2010/04/30 """
 
import matplotlib as mpl
import numpy as np
import datetime

def datenum2datetime(dnum):
    floordnum = int(np.floor(dnum)) # still a float, fromordinal needs an int
    dt = datetime.datetime.fromordinal(floordnum) + datetime.timedelta(days=dnum%1)
    return dt

def datetime2datestr(dt):
    dstr = dt.strftime('%Y/%m/%d %H:%M:%S')
    return dstr

def datenum2datestr(dnum):
    dt = datenum2datetime_alt(dnum)
    dstr = datetime2datestr(dt)
    return dstr

def datenum2datetime_alt(dnum):
    dt = mpl.dates.num2date(dnum)
    return dt

def datetime2datenum_alt(dt):
    dnum = mpl.dates.date2num(dt)
    return dnum

def epoch2datetime(et):
    """ epoch Unix times to datetime """
    epoch0 = datetime.datetime(1970,1,1,0,0,0)
    dt = epoch0 + datetime.timedelta(0, et, 0)
    return dt
        
def datetime2epoch(dt):
    """ datetime to epoch Unix time """
    epoch0 = datetime.datetime(1970,1,1,0,0,0)
    td = dt - epoch0
    et = td.days * 86400.0 + td.seconds + float(td.microseconds)/1000000
    return et
    
def utnow():
    """ epoch time now """
    dt = datetime.datetime.utcnow()
    epochnow = datetime2epoch(dt)
    return epochnow

def datetime2datenum(d):
    """ datetime to datenum """
    dnum = 366 + d.toordinal() + d.hour/24 + d.minute/3600 + d.second/86400 + d.microsecond/86400000000
    return dnum

if __name__ == "__main__":
    print("Testing utnow")
    epochnow = utnow()
    print("epoch now = %f " % epochnow)
    
    print("Testing epoch2datetime & datetime2epoch")
    dt = datetime.datetime.utcnow()
    print("datetime now is %s " % dt)
    et = datetime2epoch(dt)
    print("epoch time is %f " % et)
    dt2 = epoch2datetime(et)
    print("which is datetime %s " % dt2)
    if (dt != dt2):
        print("Problem!")
    
    print("Testing vectorized forms of datetime2epoch & epoch2datetime")
    from scipy import vectorize, r_
    vector_d2e = vectorize(datetime2epoch)
    vector_e2d = vectorize(epoch2datetime)
    dt2 = r_[[dt]*5]
    print(dt2)
    et2 = vector_d2e(dt2)
    print(et2)
