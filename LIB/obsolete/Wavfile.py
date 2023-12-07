import os
from datetime import datetime
from obspy.core import read, Stream
from obspy.io.nordic.core import blanksfile

# 9901-01-0242-12S.MVO_14_1
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
        self.filetime = datetime(int(yyyy), int(mm), int(dd), hour=int(HH), minute=int(MI), second=int(SS))

    
    def register(self, evtype, userid, overwrite=False, evtime=None):
        if not evtime:
            evtime = self.filetime
        sfilepath = blanksfile(wavfile, evtype, userid, overwrite=False, evtime=evtime)
        return sfilepath
    
    def find_sfile(self):
        """ incomplete. idea is to find the sfile, which may not have the same time """
        spath = os.path.dirname(wavfile)
        spath = spath.replace('WAV','REA',1)
        yyyy = self.evtime.year
        mm = self.evtime.month
        dd = self.evtime.day
        HH = self.evtime.hour
        MI = self.evtime.minute
        SS = self.evtime.second
        sbase = dd + "-" + HH + MI + "-" + SS + "L.S" + yyyy + mm
        sfile = os.path.join(spath, sbase)
        if not os.path.exists(sfile):
            pass        
        return sfile
               
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

