import pandas as pd
import os, glob
import obspy.core
SEISAN_DATA = os.path.join( os.getenv('HOME'),'DATA','MVO')
SEISAN_DB = 'MVOE_'

allcsv = os.path.join(SEISAN_DATA, 'reawav_%sall.csv' % SEISAN_DB)
dfall = pd.read_csv(allcsv)


# count events of each number of traces
from collections import Counter
counter_of_traces = Counter(dfall['num_traces'])
print(counter_of_traces.most_common())


# Compare compression of Seisan vs Miniseed vs Pickle files for same data
dfmaxtraces = dfall[dfall['num_traces']==36]
#print(dfmaxtraces)
print('WAV, Bzip2/Seisan, Miniseed/Seisan, Pickle/Seisan')
for i,row in dfmaxtraces.iterrows():
    WAVfile = os.path.join(SEISAN_DATA, row['path'])
    WAVfile_size = os.path.getsize(WAVfile)
    st = obspy.core.read(WAVfile)
    WAVbase = os.path.basename(WAVfile)
    MSEEDfile = WAVbase + '.mseed'
    st.write(MSEEDfile)
    MSEEDfile_size = os.path.getsize(MSEEDfile)
    PICKLEfile = WAVbase + '.pickle'
    st.write(PICKLEfile)
    PICKLEfile_size = os.path.getsize(PICKLEfile)
    os.system("cp %s %s" % (WAVfile, WAVbase))
    os.system("bzip2 %s" % WAVbase)
    BZ2file = WAVbase + ".bz2"
    BZ2file_size = os.path.getsize(BZ2file)
    print('%s, %.2f, %.2f, %.2f' % (WAVbase, BZ2file_size/WAVfile_size, MSEEDfile_size/WAVfile_size,PICKLEfile_size/WAVfile_size))
