#!/usr/bin/env python
################################
# General setup                #
################################

import os
import pandas as pd
import numpy as np
import glob


### PREPARE THE CATALOG DataFrame ###
SEISAN_DATA = os.path.join( os.getenv('HOME'),'DATA','MVO') # e.g. /home/user/seismo
pandaSeisDir = os.path.join(SEISAN_DATA, 'miniseed_c') # e.g. /home/user/seismo/pandaSeis
SEISAN_DB = 'MVOE_' # e.g. the seisan database name (e.g. MVOE_)
PROJECTDIR = os.path.join(os.getenv('HOME'),'src', 'kitchensinkGT', 'PROJECTS', 'MontserratML') # this dir
csvfile_internal = 'catalog/30_MVO_labelled_events_filtered.csv' # has to match that in AAA-master/config/general/newsettings_10.json
csvfile_internal = './AAA-master/MONTSERRAT/' + csvfile_internal

####
# Find all files matching
pattern = csvfile_internal.replace('.csv','_predicted_*.csv')
allfiles = glob.glob(pattern)
frames = []
for thisfile in allfiles:
    frames.append(pd.read_csv(thisfile))
alldf = pd.concat(frames)
lod = []


# get list of all unique corrected_DSN_mseed
unique_mseed = list(set(list(alldf['corrected_DSN_mseed'])))
print(unique_mseed)

i = 0
for this_mseed in unique_mseed:
    i+=1
    print('Processing %d of %d: %s' % (i, len(unique_mseed), this_mseed)) 
    
    thiseventdf = alldf[alldf['corrected_DSN_mseed']==this_mseed]
    eventdict = {}
    for col in thiseventdf.columns:
        for index in thiseventdf.index:
            eventdict[col] = thiseventdf.loc[index, col]
    eventdict['best_predicted_class'] = ''
    eventdict['all_predicted_classes'] = ''
    predicted_threshold = 0.5
    for classcol in ['r', 'e', 'l', 'h', 't']:
        colname = 'prob_%s' % classcol
        eventdict[classcol] = thiseventdf[colname].mean()
        if eventdict[classcol]>=predicted_threshold:
            eventdict['best_predicted_class'] = classcol
    #for row in thiseventdf:
    #    eventdict['all_predicted_classes'] = '%s%s' % (eventdict['all_predicted_classes'],row['predicted_class'])
    lod.append(eventdict)
alleventsdf = pd.DataFrame(lod)
alleventsdf.drop_duplicates(inplace=True)
alleventsdf.to_csv('30_MVO_labelled_events_filtered_predicted_networkaveraged.csv')
   