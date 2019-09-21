import os, sys, datetime
SEISAN_DATA_ROOT = "/raid/data/seisan"
DB = "MVOE_"

class Sfile:
    'Base class for Sfile parameters'
    sfilePath = ""
    year = 0
    month = 0
    day = 0
    hour = 0
    minute = 0
    second = 0.0
    startTime = datetime.datetime(1, 1, 1, 0, 0, 0)
    mainclass = ""
    subclass = ""
    op = ""
    id = ""
    sfilepath = ""
    wavfilePaths = []
    aeffilePaths = []
    
    def __init__(self, sfile):
        print "* Creating Sfile object in Python for %s *" % sfile
        file=open(sfile,'r') 
        self.sfilepath = sfile.strip()

        row = file.readlines()

        for line in row:

            if len(line)<80:
                continue

            if line[79]=='1':
                self.year = line[1:5]
                self.month = line[6:8]
                self.day = line[8:10]
                self.hour = line[11:13]
                self.minute = line[13:15]
                self.second = line[15:20]
                self.mainclass = line[21:23]
                self.startTime = datetime.datetime(int(self.year), int(self.month), int(self.day), int(self.hour), int(self.minute), int(float(self.second)))
        
            if line[79]=='I':
                self.op = line[30:33]
                self.id = line[60:74]

            if line[79]=='6':
                wavfile = line[1:36]
                wavfullpath = "%s/WAV/%s/%04d/%02d/%s" % (SEISAN_DATA_ROOT, DB, int(self.year), int(self.month), wavfile)
                self.wavfilePaths.append(wavfullpath.strip())
                aeffullpath = "%s/WAV/%s/%04d/%02d/%s" % (SEISAN_DATA_ROOT, DB, int(self.year), int(self.month), wavfile)
                if os.path.exists(aeffullpath):
                    self.aeffilePaths.append(aeffullpath.strip())
    
            if line[79]=='3':
                if line[1:10]=="VOLC MAIN":
                    self.subclass = line[11]

    def __str__(self):
        str = "S-file path: " + self.sfilepath
        str += "\n\tStart time: " + self.startTime.strftime("%Y-%m-%d %H:%M:%S")
        str += "\n\tClass: " + self.mainclass
        str += "\n\tSubclass: " + self.subclass
        str += "\n\tAnalyst: " + self.op
        str += "\n\tId: " + self.id
        for wavfullpath in self.wavfilePaths:
            if os.path.exists(wavfullpath):
                str += "\n\twavfile: " + wavfullpath + " found"
            else:
                str += "\n\twavfile: " + wavfullpath + " missing"                    
        for aeffullpath in self.aeffilePaths:
            if os.path.exists(aeffullpath):
                str += "\n\taeffile: " + aeffullpath

    def export(self):

class Wavfile(self):
    def __init__(self, path=None):
        self.path = path.strip()
    self.st = None

    def plot(self):
       #if os.path.exists(self.path):
       #     sys.argv = [mulplt.__file__, self.path]
       #     mulplt.main()
       if self.st:
           self.st.plot()
       else:
           self.read()
           self.st.plot()

    def mulplt(self):
       if self.st:
           self.st.mulplt()
       else:
           self.read()
           self.st.mulplt()

    def read(self):
        #st = None
        self.st = StreamGT()
        if os.path.exists(self.path):
            self.st.read(self.path)

    def __str__(self):
        str = "wavfile: %s %r " % (self.path, os.path.exists(self.path))
        if self.st:
            str += "\n\t" + self.st.__str__()

class AEFfile(self):
    def __init__(self, path=None):
        print "* Creating AEFfile object in Python for %s *" % self.path
        self.path = path.strip()
        lines = list()
        if os.path.exists(path):
            file=open(path,'r') 
            lines = file.readlines()
        for line in lines:

            if len(line)<80:
                continue

            if line[79]=='1':
                self.year = line[1:5]
                self.month = line[6:8]
                self.day = line[8:10]
                self.hour = line[11:13]
                self.minute = line[13:15]
                self.second = line[15:20]
                self.mainclass = line[21:23]
                self.startTime = datetime.datetime(int(self.year), int(self.month), int(self.day), int(self.hour), int(self.minute), int(float(self.second)))
        
            if line[79]=='I':
                self.op = line[30:33]
                self.id = line[60:74]

            if line[79]=='6':
                wavfile = line[1:36]
                wavfullpath = "%s/WAV/%s/%04d/%02d/%s" % (SEISAN_DATA_ROOT, DB, int(self.year), int(self.month), wavfile)
                self.wavfilePaths.append(wavfullpath.strip())

    def __str__(self):
        str = "aeffile: %s %r " % (self.path, os.path.exists(self.path))
            
if __name__ == "__main__":

    # read in a single Sfile
    sfilepath = ''
    s = Sfile(sfilepath)

    # inspect it
    print s

    # export to ObsPy event object
    e = s.export('obspy')
    print e

