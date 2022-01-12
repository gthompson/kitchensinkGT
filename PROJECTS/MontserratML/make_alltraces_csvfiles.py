import os
from glob import glob
SEISAN_DATA = '/raid/data/Montserrat/MASTERING/seisan'
SEISAN_DB = 'MVOE_'
# Make the monthly all_traces PD files    
DBDIR = os.path.join(SEISAN_DATA, 'miniseed_c',SEISAN_DB)
os.chdir(DBDIR)
yeardirs = sorted(glob('[12]???'))
for yeardir in yeardirs:
    print(yeardir)
    YYYY = os.path.basename(yeardir)
    monthsdirs = sorted(glob(os.path.join(yeardir,'[01]?')))
    for monthdir in monthsdirs:
        MM = os.path.basename(monthdir)
        print(monthdir)
        os.system('cat %s/%s/*.csv | sort | uniq > alltraces_%s%s.csv' % (YYYY, MM, YYYY, MM) )
        

# Make the full all_traces PD file
os.system('cat alltraces_??????.csv | sort | uniq > alltraces.csv')
