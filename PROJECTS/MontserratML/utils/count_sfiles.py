#!/usr/bin/env python
import glob
import os
from collections import Counter
total_DSN_traces = 0
total_ASN_traces = 0
all_sfiles = []
all_DSN_wavfiles = []
all_ASN_wavfiles = []
all_subclasses = []
HOME = os.getenv('HOME')
os.chdir(os.path.join(HOME, 'DATA', 'MVO'))
print('Directory, cum_N_sfiles, cum_N_DSN_wavfiles, cum_N_DSN_traces, cum_N_ASN_wavfiles, cum_N_ASN_traces')
for ydir in glob.glob('REA/MVOE_/????'):
    YYYY = ydir[-4:]
    for mdir in glob.glob('%s/??' % ydir):

        sfiles = glob.glob('%s/*.S??????' % mdir)
        MM = mdir[-2:]

        for this_sfile in sfiles:
            all_sfiles.append(this_sfile)
            fptr=open(this_sfile, 'r')
            Lines = fptr.readlines()
            subclass = None
            if 'R.S' in this_sfile:
                subclass = 'REGIONAL'
            if 'D.S' in this_sfile:
                subclass = 'TELESEISM'
            if 'L.S' in this_sfile:
                subclass = 'LOCAL'

            # Strips the newline character
            for line in Lines:
                line = line.strip()
                if len(line)<75:
                    continue
                if line[-1] == '3' and 'VOLC MAIN' in line:
                    subclass = line[10]
                if line[-1] == '6':
                    this_wavfile = line[0:70].strip()
                    if this_wavfile[0]=='M':
                        continue
                    num_traces = 0
                    if this_wavfile[-2]=='_':
                        num_traces = int(this_wavfile[-4:-2])
                    elif this_wavfile[-5:-3] == '__':
                        num_traces = int(this_wavfile[-2:])
                    #print(this_wavfile, num_traces)
                    if 'MVO' in this_wavfile:
                       total_DSN_traces += num_traces
                       all_DSN_wavfiles.append(this_wavfile)
                    if 'SPN' in this_wavfile:
                       total_ASN_traces += num_traces 
                       all_ASN_wavfiles.append(this_wavfile)
            fptr.close()
            all_subclasses.append(subclass)
        print('%s, %7d, %7d, %7d, %7d, %7d' %(mdir, len(all_sfiles), len(all_DSN_wavfiles), total_DSN_traces, len(all_ASN_wavfiles), total_ASN_traces))
        counter_of_classes = Counter(all_subclasses)
        #print(YYYY, MM, counter_of_classes.most_common())
#counter_of_classes = Counter(all_subclasses)
#print(counter_of_classes.most_common())
