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
