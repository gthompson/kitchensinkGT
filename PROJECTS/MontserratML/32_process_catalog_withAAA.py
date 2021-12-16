#!/usr/bin/env python
################################
# General setup                #
################################

import os
import pandas as pd
import pickle
import numpy as np


### PREPARE THE CATALOG DataFrame ###
SEISAN_DATA = os.path.join( os.getenv('HOME'),'DATA','MVO') # e.g. /home/user/seismo
pandaSeisDir = os.path.join(SEISAN_DATA, 'miniseed_c') # e.g. /home/user/seismo/pandaSeis
SEISAN_DB = 'MVOE_' # e.g. the seisan database name (e.g. MVOE_)
PROJECTDIR = os.path.join(os.getenv('HOME'),'src', 'kitchensinkGT', 'PROJECTS', 'MontserratML') # this dir
#csvfile_external = os.path.join(PROJECTDIR, 'MVO_labelled_events.csv')
csvfile_external = os.path.join(SEISAN_DATA, 'MachineLearning', SEISAN_DB, 'runAAA', 'MVOE_11_labelled_events.csv')
#csvfile_internal = './catalog/MVO_labelled_events_filtered.csv'
csvfile_internal = 'catalog/30_MVO_labelled_events_filtered.csv' # has to match that in AAA-master/config/general/newsettings_10.json
csvfile_internal = './AAA-master/MONTSERRAT/' + csvfile_internal
output_path_cat = csvfile_internal.replace('.csv', '.pd')
alltraces_file = '30_alltraceDFs.csv'

################################
# Machine learning and testing #
################################
import sys
sys.path.insert(0, './AAA-master/automatic_processing')
#import tools
from config import Config
from analyzer import Analyzer

# Change if you want your screen to keep quiet
# 0 = quiet
# 1 = in between
# 2 = detailed information
verbatim = 1

# Init project with configuration file
config = Config('./AAA-master/config/general/newsettings_10.json', verbatim=verbatim)
config.readAndCheck()  

##########################
# Variables to loop over #
##########################
alltraces = pd.read_csv(alltraces_file)

traceIDs = ['MV.MBWH..SHZ', 'MV.MBLG..SHZ', 'MV.MBRY..SHZ']
minWeights = range(4)
classes_to_include = [ ['l', 't'], ['e', 'r'], ['h', 'l', 't'], ['h', 'l', 't', 'r'], ['e', 'h', 'l', 't', 'r'] ]

traceIDs = ['MV.MBWH..SHZ', 'MV.MBLG..SHZ', 'MV.MBGB..BHZ', 'MV.MBGH..BHZ']
minWeights = range(4)
classes_to_include = [ ['r', 'e', 'l', 'h', 't']]
combine_re = False

traceIDs = ['MV.MBWH..SHZ', 'MV.MBGB..BHN', 'MV.MBGB..BHZ', 'MV.MBGB..BHE', 'MV.MBLG..SHZ', 'MV.MBGH..BHZ', 'MV.MBGH..BHE', 'MV.MBGH..BHN', 'MV.MBGA..BHE', 'MV.MBGA..BHZ', 'MV.MBGA..BHN', 'MV.MBGE..BHE', 'MV.MBGE..BHZ', 'MV.MBGE..BHN', 'MV.MBBE..BHE', 'MV.MBBE..BHZ', 'MV.MBBE..BHN']
minWeights = [3]
classes_to_include = [ ['r', 'e', 'l', 'h', 't'] ]
combine_re = False

#############
# Functions #
#############

def cat_filter_traceID(cat, alltraces, traceID):
    # subset catalog based on traceID
    matchingEvents = alltraces[alltraces['id']==traceID]
    for i, row in cat.iterrows():
        matching_indices = matchingEvents.index[matchingEvents['filetime']==row['filetime']].tolist()
        if len(matching_indices)==1:
            pass
        else:
            cat.drop(i, inplace=True)
    N = len(cat.index)
    #print('%d events after matching against traceID' % N)
    return N

def cat_filter_classes(cat, remove_classes):
    """
    for rmclass in remove_classes:
        print('Removing %s' % rmclass)
        cat = cat[cat['class']!=rmclass]
    """
    cat=cat[cat["class"].isin(remove_classes)]
    N = len(cat.index)
    #print('%d events after removing classes' % N)
    return cat, N

def cat_filter_weight(cat, minWeight):
    #if minWeight>0:
    #    cat = cat[cat['weight']>=minWeight]
    cat=cat[cat["weight"].isin(range(minWeight,13))]
    N = len(cat.index)
    #print('%d events after filtering above %d' % (N, minWeight))
    return cat, N

def cat_check_numbers(cat, minthresh = 20):
    df = cat.copy()
    tooSmall = False
    lengths = []
    for subclass in df['class'].unique():
        dfs = df[df['class']==subclass]
        N = len(dfs.index)
        if N<minthresh:
            tooSmall=True
        lengths.append(N)
    return tooSmall, lengths


#######################
# Looping starts here #
#######################
minPerClass = 30
counter = 0


for combine_re in [False, True]:
    results_list = []

    for traceID in traceIDs:

        # what traceID are we looking for - read_montserrat needs this - should figure out how to write this into the config
        fptr = open('./AAA-master/MONTSERRAT/current_traceID.txt','w')
        fptr.write(traceID)
        fptr.close()

        for minWeight in minWeights:
            for include_classes in classes_to_include:
                # reload cat because we filter it down each time
                #cat = pickle.load(open(output_path_cat,'rb'))
                cat = pd.read_csv(csvfile_internal)
                if combine_re:
                    cat.loc[cat['class']=='e', 'class']='r'
                    if 'e' in include_classes:
                        include_classes.remove('e')

                print(cat['class'].value_counts())
                print(traceID, include_classes, minWeight)

                results_dict = {}
                results_dict['traceID'] = traceID
                results_dict['classes'] = ','.join(include_classes)
                results_dict['minWeight'] = minWeight
                results_dict['NtraceID'] = cat_filter_traceID(cat, alltraces, traceID)
                cat, N = cat_filter_classes(cat, include_classes) 
                results_dict['Nclasses'] = N
                cat, N = cat_filter_weight(cat, minWeight)
                results_dict['Nweight'] = N
                tooSmall, lengths = cat_check_numbers(cat, minPerClass)
                results_dict['counts'] = str(lengths)

                results_dict['acc_mean'] = None
                results_dict['acc_std'] = None

                if N>=minPerClass*len(include_classes) and not tooSmall:
                    #try:
                        print(cat.groupby('class').size())
                        analyzer = Analyzer(config, verbatim=verbatim, catalog=cat)
                        allData, allLabels, acc, allPredictions, allProbabilities = analyzer.learn(config, returnData=True, verbatim=verbatim) # If you want the data
                        results_dict['acc_mean'] = np.round(np.mean(acc)*100, 1)
                        results_dict['acc_std'] = np.round(np.std(acc)*100, 1)                        
                        cat['predicted_class'] = allLabels
                        cat['traceID'] = traceID
                        for i, classcol in enumerate(sorted(include_classes)):
                            colname = 'prob_%s' % classcol
                            cat[colname] = allProbabilities[:,i]
                        cat.to_csv(csvfile_internal.replace('.csv','_predicted_%d.csv' % counter)) 
                results_list.append(results_dict)
                counter += 1

    resultsDF = pd.DataFrame(results_list)
    resultsCSV = ''.join(include_classes) + '.csv'
    resultsDF.to_csv(resultsCSV)
