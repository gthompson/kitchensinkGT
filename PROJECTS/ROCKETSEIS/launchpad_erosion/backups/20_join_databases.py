#!/usr/bin/env python
# coding: utf-8

get_ipython().run_line_magic('run', 'header.ipynb')
#print(paths)

os.chdir(paths['outdir'])
import datetime
from IPython.display import clear_output
#ANTELOPE = os.getenv('ANTELOPE')
years = ['2022']
#nets = ['6I', '6S', 'FL', 'AM', 'XA']
nets = ['FL', 'AM', 'XA']
YYYY = '2022'

def dbcreate(dbname, schema="css3.0", dbpath=None):
    nl = '\n'
    dbpaths = ""
    if dbpath:
        for p in dbpath:
            d = os.path.dirname(p)
            b = os.path.basename(p)
            if dbpaths:
                connector = ":"
            else:
                connector = ""
            #thisdbpath = "%s/{%s}" % (d,b) # no, in same directory as descriptor, so just need dbbasename
            thisdbpath = "{"+b+"}"
            dbpaths += connector + thisdbpath
    contents = "#\n"
    contents += f"schema {schema}{nl}"
    if dbpath:
        contents += f"dbpath {dbpaths}{nl}"
    print(contents)
    with open(dbname, "w") as f:
        f.write(contents)

for YYYY in years:
    for jday in range(83,367): # 83 is first day of seismic data in SDS 
        start_date = datetime.datetime.strptime(f"{YYYY}{jday:03d}", '%Y%j')
        print(start_date, end="\n")
        ymd = start_date.strftime("%Y%m%d")
        #startepoch = int(start_date.strftime("%s"))
        #endepoch = starttime + 86400
        dbwellday = f"db/dbgood{ymd}"
        dbseismo = f"db/dbseismo{ymd}"
        for net in nets:
            mseed_dirs = sorted(glob.glob(os.path.join('SDS', YYYY, net, '???*', '[HD]??.D')))
            for mseed_dir in mseed_dirs:
                mseedfiles = sorted(glob.glob(os.path.join(mseed_dir, '%s.*.D.*.%03d' % (net, jday) ) ))    
                if len(mseedfiles)>0:
                    mseedfilesstr = " ".join(mseedfiles)
                    #os.system("bash ~/Developer/miniseed2days_wrapper.sh %s miniseed/good" % start_date)
                    os.system(f"miniseed2db {mseedfilesstr} {dbseismo} ")
        # create descriptor
        if os.path.isfile(dbwellday+'.wfdisc') or os.path.isfile(dbseismo+'.wfdisc'):
            dball = f"db/dberosion{ymd}"
            dbcreate(dball, schema='css3.0', dbpath=[dbwellday, dbseismo])
        clear_output(wait=True)


