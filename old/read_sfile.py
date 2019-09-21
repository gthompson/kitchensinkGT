#!/raid/apps/OBSPY/bin/python
# Read S-file (from REA)
# Glenn Thompson 2013/01/06
SEISAN_DATA_ROOT = "/raid/data/seisan"
DB = "MVOE_"
import os, sys, datetime

# read in the wav file name from the command line
sfile = sys.argv[1]

file=open(sfile,'r') 

row = file.readlines()
wavfile = []
for line in row:
    #print line[79]
    if line[79]=='1':
        year = line[1:5]
        month = line[6:8]
        day = line[8:10]
        hour = line[11:13]
        minute = line[13:15]
        second = line[15:20]
        mainclass = line[21:23]
        start_time = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(float(second)))
        print(start_time.strftime("%Y-%m-%d %H:%M:%S"), mainclass)
    if line[79]=='I':
        op = line[30:33]
        sfile_id = line[60:74]
        print(op, sfile_id)
    if line[79]=='6':
        wavfile.append(line[1:36])
        print(wavfile[-1])
        wavfullpath = "%s/WAV/%s/%04d/%02d/%s" % (SEISAN_DATA_ROOT, DB, int(year), int(month), wavfile[-1])
        print(wavfullpath)
        if os.path.exists(wavfullpath):
            print(wavfullpath + " found")
        else:
            print(wavfullpath + " missing") 
    if line[79]=='3':
        if line[1:10]=="VOLC MAIN":
            subclass = line[11]
            print(subclass) 

        

        

