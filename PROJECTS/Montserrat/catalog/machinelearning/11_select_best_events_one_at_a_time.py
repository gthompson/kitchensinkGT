#!/usr/bin/env python
# coding: utf-8

# # Montserrat event selector for Machine Learning
# The aim of this code is to find the best N events of each type, and create a corresponding CSV file and data structure for entry into Alexis' and Marielle's AAA codes.

##########################################################
# IMPORTS BEGIN ##########################################

import os
from glob import glob
import pandas as pd
import numpy as np
import sys
LIBpath = os.path.join( os.getenv('HOME'),'src','kitchensinkGT', 'LIB')
sys.path.append(LIBpath)
from libseisGT import add_to_trace_history, plot_seismograms 
from modutils import yn_choice
#from IPython.display import clear_output
from obspy import read, read_inventory #, remove_response
from libMVO import fix_trace_id, inventory_fix_id_mvo, load_mvo_inventory, plot_station_amplitude_map, read_enhanced_stream
from metrics import check_for_spikes
#from statsmodels.stats.weightstats import DescrStatsW
#from obspy.core import read
#import matplotlib.pyplot as plt
CWD = os.getcwd() # might need to pass explicitly to qc_event()
sys.path.append(CWD)

# IMPORTS END ############################################
##########################################################
# FUNCTIONS BEGIN ########################################

def read_volcano_def(volcanodefcsv):
    subclass_df = pd.read_csv(volcanodefcsv)
    subclass_df.columns = subclass_df.columns.str.lstrip()
    return subclass_df

#

#def build_master_event_catalog(csvdir, seisandbname, catalogfile, subclasses_for_ML, max_duration = 300):
def build_master_event_catalog(pandaSeisDBDir, SEISAN_DB, catalogfile, subclasses_for_ML, max_duration = 300):
    # load all the year/month CSV files
    #csvfiles = glob(os.path.join(pandaSeisDBDir, 'reawav_%s??????.csv' % SEISANDB))
    #csvfiles = glob(os.path.join(pandaSeisDBDir, 'reawav_%s[12][0-9][0-9][0-9][0-1][0-9].csv' % SEISAN_DB))
    csvfiles = glob(os.path.join(pandaSeisDBDir, 'catalog_%s[12][0-9][0-9][0-9][0-1][0-9].csv' % SEISAN_DB))
    frames = []
    if len(csvfiles)==0:
        print('No reawav*.csv files found. Cannot proceed')
        exit()
    for csvfile in csvfiles:
        df = pd.read_csv(csvfile)
        frames.append(df) 
    dfall = pd.concat(frames, sort=True)
    # perhaps path should be the index instead? no, just let it go through without indexing or sorting
    #dfall.set_index('filetime', inplace=True) # we will need this later to remerge
    #dfall.set_index('path', inplace=True)
    #dfall.sort_index(inplace=True)
    """
    for index, row in dfall.iterrows():
        
        # For simplicity, copy 'D' and 'R' mainclass to subclass
        if row['mainclass'] in ['D', 'R']: # Do I need L here too?
            dfall.loc[index, 'subclass']=row['mainclass']
    """ 
    # replace loop above
    for mainclass in ['R', 'D']:
        dfall.loc[dfall['mainclass'] == mainclass, 'subclass'] = mainclass
    
    # Drop the mainclass column, as it is now superfluous.
    dfall.drop(columns=['mainclass'], inplace=True)
    
    # Add an etype column
    #dfall['etype'] = dfall['subclass'].replace(subclasses, etypes)
    
    # Add columns to assign a percentage for each subclass
    for subclass in subclasses_for_ML:
        dfall[subclass] = 0
    
    # But set column for actual subclass to 100%  
    """
    for index, row in dfall.iterrows():            
        # Set dfall['h'] to 100 if dfall['subclass']=='h', etc
        dfall.loc[index, row['subclass']] = 100        
    """
    # replace row operations above
    for subclass in subclasses_for_ML:
        dfall.loc[dfall['subclass'] == subclass, subclass] = 100
        
    # Add a new_subclass column
    dfall['new_subclass'] = dfall['subclass']

    # Add weight column. I will give really clear events higher weight when I process them
    dfall['weight']=3 # weight for events I have not visually checked
    
    # Add column that records if event is checked
    dfall['checked']=False
    
    # Add column that records if event is marked for splitting
    dfall['split']=False    
    
    # Add column that records if event is marked for deletion
    dfall['delete']=False
    
    # Add column that records if event should be ignored
    # Ignore any events longer than 1-minute, as they are likely to contain multiple events 
    # or just be unhelpful for classifying short signals which are more common
    # SCAFFOLD - the twin column no longer seems to exist
    #dfall['ignore'] = dfall['twin']>max_duration
    dfall['ignore'] = False
    
    # Now we have a catalog dataframe we can work with. Let's save this.
    #dfall2 = dfall.reset_index()    
    #dfall2.to_csv(catalogfile)
    _df2file_without_index(dfall.copy(), catalogfile)
    
    return dfall

#    
    
def _select_best_events(df, allowed_subclasses, N=100, exclude_checked=True):
    # When we are iterating, we want to exclude the checked events, else we are repeating our work.
    # When we return our final list, we do not want to exclude checked events
    all_subclasses = df['new_subclass'].unique()
    best_events_dict = {}
    for subclass in all_subclasses:
        if subclass in allowed_subclasses:
            #print('\nProcessing subclass=%s' % subclass)
            is_subclass =  df['new_subclass']==subclass
            df_subclass = df[is_subclass]
            
            # mechanism to weight out checked events
            df_subclass['include'] = 1 - df_subclass['ignore']
            if exclude_checked:
                print('Excluding checked events')
                df_subclass['include'] = df_subclass['include'] * (1 - df_subclass['checked'])
            
            if len(df_subclass.index)>0:
                # we use three criteria for ranking events. detection_quality has the largest magnitude, but can be missing, so we also add snr, 
                # and finally quality as a tie-braker, since it has a small range for events that have made it this far
                df_subclass['sortcol'] = df_subclass['quality'] # always present
                if 'snr' in df_subclass.columns:
                    df_subclass['sortcol'] = df_subclass['sortcol'] + df_subclass['snr']
                if 'detection_quality' in df_subclass.columns:
                    df_subclass['sortcol'] = df_subclass['sortcol'] + df_subclass['detection_quality']                    
                df_subclass['sortcol'] = df_subclass['sortcol'] * df_subclass['weight'] * df_subclass['include']
                df_subclass.drop(columns=['include'], inplace=True)

                L = len(df_subclass.index)
                H = int(min([L, N]))
                print('Selecting %d events of type %s from a total of %d' % (H, subclass, L))
                df_subclass.sort_values(by=['sortcol'], ascending=False, inplace=True)
                df_subclass.drop(columns=['sortcol'], inplace=True)
                df_subclass = df_subclass.head(H)
                #if 'bandratio_[1.0_6.0_11.0]' in df_subclass.columns:
                #    df_subclass.rename(columns = {'bandratio_[1.0_6.0_11.0]':'band_ratio'}, inplace = True)
                best_events_dict[subclass]=df_subclass

    return best_events_dict

#

def get_fingerprints(dfall, allowed_subclasses, N=100, exclude_checked=False):
    """
    All we do right now is a dataframe describe, so we return the stats of each column.
    """
    #df = best_events_dict.groupby("subclass")
    fingerprints = {}
    best_events_dict = _select_best_events(dfall, allowed_subclasses, N=N, exclude_checked=exclude_checked)
    for subclass in allowed_subclasses:
        if subclass in best_events_dict.keys(): 
            #print('Computing fingerprint for subclass ',subclass)
            #df_subclass = df.get_group(subclass)
            df_subclass = best_events_dict[subclass]     
            fingerprints[subclass] = df_subclass[[ 'peaktime', 
                'kurtosis', 'medianF', 'peakF', 'bw_min', 'bw_max', 'bandratio_[1.0_6.0_11.0]']].describe()
            #print(fingerprints[subclass])
        else:
            if os.path.exists('fingerprint_%s.csv' % subclass):
                fingerprints[subclass]=pd.read_csv('fingerprint_%s.csv' % subclass)
    return fingerprints

"""    
def get_weighted_fingerprints(dfall, subclasses_for_ML, N=300, exclude_checked=False):

    fingerprints = {}
    best_events_dict = _select_best_events(dfall, subclasses_for_ML, N=N, exclude_checked=exclude_checked)
    for subclass in subclasses_for_ML:
        if subclass in best_events_dict.keys(): 
            print('Computing fingerprint for subclass ',subclass)
            thisdf = best_events_dict[subclass]
            if len(thisdf.index)>30:
                statsdf = pd.DataFrame()
                statsdf['statistic'] = ['mean', 'std', '25%', '50%', '75%']
                statsdf.set_index(['statistic'], inplace = True)

                for col in [ 'peaktime', 'kurtosis', 'medianF', 'peakF', 'bw_min', 'bw_max', 'band_ratio']:
                    # compute mean, std, median, 25% percentile, 75% percentile
                    wdf = DescrStatsW(thisdf[col], weights=thisdf['weight'].astype(float)*thisdf[subclass].astype(float), ddof=1)
                    p = [0.25,0.50,0.75]
                    q  = wdf.quantile(p) 
                    statsdf.loc['mean', col] = wdf.mean
                    statsdf.loc['std', col] = wdf.std
                    statsdf.loc['50%', col] = q[p[1]]
                    statsdf.loc['25%', col] = q[p[0]]
                    statsdf.loc['75%', col] = q[p[2]]

                fingerprints[subclass] = statsdf
                print(fingerprints[subclass])
    return fingerprints     
"""
 
#    

def save_fingerprints(fingerprints, allowed_subclasses):
    for subclass in allowed_subclasses:
        if subclass in fingerprints.keys(): 
            _df2file_without_index(fingerprints[subclass], 'fingerprint_%s.csv' % subclass)
            
#

def select_next_event(dfall, subclasses_for_ML):
    """ The goal here is just to pick the next event of a particular subclass from the unchecked events """
    checked = dfall[dfall['checked']==True]
    unchecked = dfall[dfall['checked']==False]
    subclasses = subclasses_for_ML.copy()
    
    tryagain = True
    while tryagain and len(subclasses)>0:
    
        # Check how many we have of each class
        counts = None
        nextclass = 'r'
        mincounts = 99999

        for subclass in subclasses:
            dfs = checked[checked['new_subclass']==subclass]
            counts=len(dfs.index)
            if counts<mincounts:
                mincounts=counts
                nextclass=subclass
        #print('next event is of subclass = %s' % nextclass)

        is_subclass = unchecked['subclass']==nextclass
        df = unchecked[is_subclass]

        # mechanism to weight out checked events
        df = df[df['ignore']==False]
        
        # weight out events with less than 3 traces
        df = df[df['num_traces']>2]

        if len(df.index)>0:
            # we use three criteria for ranking events. detection_quality has the largest magnitude, but can be missing, so we also add snr, 
            # and finally quality as a tie-braker, since it has a small range for events that have made it this far
            df['sortcol'] = df['quality'] # always present
            if 'snr' in df.columns:
                df['sortcol'] = df['sortcol'] + df['snr']
            if 'detection_quality' in df.columns:
                df['sortcol'] = df['sortcol'] + df['detection_quality']                    
            df['sortcol'] = df['sortcol'] * df['weight']
            #df.drop(columns=['include'], inplace=True)
            df.sort_values(by=['sortcol'], ascending=False, inplace=True)
            df.drop(columns=['sortcol'], inplace=True)
            #df = df.head(1)
            #if 'bandratio_[1.0_6.0_11.0]' in df.columns:
            #    df.rename(columns = {'bandratio_[1.0_6.0_11.0]':'band_ratio'}, inplace = True)
            return df.index[0] # a WAVfile path
        else:
            subclasses.remove(nextclass)
    return None

#
"""
from obspy.core import read
def _deconvolve_instrument_response(st):
    for tr in st:
        this_inv = None
        need_to_correct = False
        if 'units' in tr.stats:
            if 'units'=='Counts':
                need_to_correct = True
        else:
            need_to_correct = True
        if need_to_correct:           
            if bool_correct_data: # try to find corresponding station XML
                this_inv = load_mvo_inventory(tr, os.path.join(SEISAN_DATA, 'CAL'))
        if this_inv and tr.stats['units'] == 'Counts':
            tr.remove_response(inventory=this_inv, output="VEL")
            tr.stats['units'] = 'm/s'
            add_to_trace_history(tr, 'deconvolved')
"""
#

def _guess_subclass(row, fingerprints, subclasses_for_ML):
    chance_of = {}
    #print(row)
    for item in subclasses_for_ML:
        chance_of[item]=0.0
    for subclass in fingerprints.keys():
        fingerprint = fingerprints[subclass]
        for param in fingerprint.columns:
            thisval = row[param] # the value of this parameter for the event
            
            # test against mean+/-std
            meanval = fingerprint.loc['mean', param]
            stdval = fingerprint.loc['std', param]
            minus1sigma = meanval - stdval
            plus1sigma = meanval + stdval
            
            if thisval > minus1sigma and thisval < plus1sigma:
                weight = 1.0 - abs(thisval-meanval)/stdval
                chance_of[subclass] += weight
                
            # test against 25-75% percentile
            medianval = fingerprint.loc['50%',param]
            val25 = fingerprint.loc['25%',param]
            val75 = fingerprint.loc['75%',param]
            if thisval > val25 and thisval < val75:
                if thisval < medianval:
                    weight = 1.0 - (medianval-thisval)/(medianval-val25)
                else:
                    weight = 1.0 - (thisval-medianval)/(val75-medianval)
                chance_of[subclass] += weight
    
    print('The event is classified as %s, but here are our guesses:' % row['subclass'])
    total = 0
    for subclass in chance_of.keys():
        total += chance_of[subclass]
    for subclass in chance_of.keys():
        if total>0:
            print('%s: %.0f' % (subclass, 100*chance_of[subclass]/total), end=', ')


from pathlib import Path
def qc_event(df, thiswav, subclasses_for_ML, seisan_subclasses, fingerprints, pandaSeisDBDir, station0hypfile):
    
    subclass = df.loc[thiswav, 'new_subclass']
    orig_subclass = df.loc[thiswav, 'subclass']
    mseedpath = df.loc[thiswav, 'corrected_DSN_mseed']
    empty_df = pd.DataFrame()
    p = Path(mseedpath)
    mseedpath = os.path.join(pandaSeisDBDir, p.parts[-3], p.parts[-2], p.parts[-1])

    
    # we are adding a check here for spikes that will now run as part of 00_ but wasn't in place when we ran it before
    # so we add it here. if we suspect spikes, we mark that trace ID for removal and do this after we load the MSEED file below. 
    # The reason we do the spike check on the raw WAV file is because filters run to produce the MSEED file distort the spike
    if os.path.exists(mseedpath):
        wavpath = mseedpath.replace('miniseed_c', 'WAV').replace('.mseed','')
        rawst = read(wavpath)
        for tr in rawst:
            check_for_spikes(tr)
            
        good_traces = 0
        trace_ids_to_eliminate = []
        fix_trace_id(rawst)
        for tr in rawst:
            check_for_spikes(tr)
            if tr.stats.quality_factor > 1.0:
                good_traces += 1
            else:
                trace_ids_to_eliminate.append(tr.id)


        # plot whole event file
        #st = read(mseedpath) #.select(station='MBWH')
        st = read_enhanced_stream(mseedpath)
        
        # adding due to spike checks above
        for tr in st:
            if tr.id in trace_ids_to_eliminate:
                st.remove(tr)
        if len(st)<3:
            print('not enough traces after removing spiky traces')
            df.loc[thiswav, 'ignore'] = True
            df.loc[thiswav, 'checked'] = True
            return False
                
        st.filter('bandpass', freqmin=0.5, freqmax=25.0, corners=4)
        #_deconvolve_instrument_response(st)

        # compute ampeng to show later
        """
        traceids = []
        energies = []
        amplitudes = []
        for tr in st:
            amp = np.max(np.absolute(tr.data))
            energy = np.sum(np.square(tr.data))/tr.stats.sampling_rate
            traceids.append(tr.id)
            amplitudes.append(amp)
            energies.append(energy)
        dfenergy = pd.DataFrame()
        dfenergy['traceID']=traceids
        dfenergy['amplitude']=amplitudes
        dfenergy['energy']=energies
        dfenergy.sort_values(by=['amplitude'],ascending=False,inplace=True)
        suggested_weight = None
        if st[0].stats.units == 'Counts':
            suggested_weight = int(np.log10(dfenergy.iloc[0]['energy']/2)) # works with uncorrected data in Counts only 
        """

        # plot all
        st.plot(equal_scale=False, outfile=os.path.join(CWD, 'current_event.png'));
        #st.plot(equal_scale=False, outfile=os.path.join(CWD, 'current_event.png'), method='full');
        #plot_seismograms(st, outfile=os.path.join(CWD, 'current_event.png'))
        # also plot a fixed 30-seconds around the peaktime
        #starttime = st[0].stats.starttime + df.iloc[0,'peaktime']-10
        starttime = st[0].stats.starttime + df.loc[thiswav,'peaktime']-10
        endtime = starttime + 20
        st.trim(starttime=starttime, endtime=endtime, pad=True, fill_value=None)
        st.plot(equal_scale=False, outfile=os.path.join(CWD, 'current_event_zoomed.png'));
        #st.plot(equal_scale=False, outfile=os.path.join(CWD, 'current_event_zoomed.png'), method='full');
        #plot_seismograms(st, outfile=os.path.join(CWD, 'current_event_zoomed.png'))

        # Show a map of amplitude distribution
        try:
            plot_station_amplitude_map(st, station0hypfile=station0hypfile, outfile='current_map.png')
        except:
            print('plot_station_amplitude_map has crashed')

        # show webpage
        if sys.platform=='linux':
            os.system("xdg-open current_event.html")
        #elif sys.platform=='darwin': # does not work, so just manually open current_event.html instead
        #    os.system("open current_event.html")

        # TEXT OUTPUT STARTS HERE
        print(' ')
        #print('Checked events currently:')
        #print(dfall[dfall['checked']==True].groupby('new_subclass').size())
        #print(row)
        #print(' ')
        print('Loaded %s, original class: %s, current class: %s' % (mseedpath, orig_subclass, subclass))

        # station amplitudes and energies
        #print(dfenergy)

        # trace df metrics
        csvpath = mseedpath.replace('.mseed', '.csv')
        tracedf = pd.read_csv(csvpath)
        #if 'bandratio_[1.0_6.0_11.0]' in tracedf.columns:
        #    tracedf.rename(columns = {'bandratio_[1.0_6.0_11.0]':'band_ratio'}, inplace = True)
        tracedf.sort_values(by=['peakamp'], ascending=False, inplace=True)
        print(' ')
        #print(tracedf[['id', 'medianF', 'bw_min', 'peakF', 'bw_max', 'band_ratio', 'kurtosis']])
        print(tracedf[['id', 'medianF', 'bw_min', 'peakF', 'bw_max', 'bandratio_[1.0_6.0_11.0]', 'kurtosis']])

        # fingerprint guesses
        print(' ')
        _guess_subclass(df.loc[thiswav], fingerprints, subclasses_for_ML)
        #print(suggested_weight)

        # Input
        checked = False
        print(' ')
        print('RECLASSIFY')
        print('Valid subclasses according to VOLCANO.DEF are: ', seisan_subclasses, end = ' ' )
        print('e.g. l, 75, h, 25, 6')
        #print('To enter percentage probabilities, e.g. 75% l, 25%h, enter l, 75, h, 25')
        #print('Optionally add a weight [0-9] too with a trailing integer, e.g. l, 75,  h, 25, 5')
        #print('Or:\n\ts = mark event for splitting')
        #print('\td = mark event for deletion')
        #print('\ti = ignore event')
        #print('\tq = quit')
        print('Or s=split, d=delete, i=ignore, q=quit')

        try:                       
            new_subclass = input('\t ?') 
            if not new_subclass:
                new_subclass = subclass
            if new_subclass == 'q': # GOOD, QUITTING
                #return empty_df, True
                return True
            if new_subclass == 's':
                  df.loc[thiswav, 'split'] = True
                  checked = True
            if new_subclass == 'i':
                  df.loc[thiswav, 'ignore'] = True
                  checked = True
            if new_subclass == 'd':
                  df.loc[thiswav, 'delete'] = True
                  checked = True
            if not checked:                         
                if not ',' in new_subclass: # convert to a subclass, percentage string
                    new_subclass = new_subclass + ', 100'
                spl = new_subclass.split(',') # split string to subclass probability list 
                if len(spl) % 2 == 1:
                    df.loc[thiswav, 'weight'] = int(spl.pop())
                spd = {spl[a].strip():spl[a + 1] for a in range(0, len(spl), 2)} # subclass probability dict
                for key in subclasses_for_ML:
                    if key in spd.keys():
                        df.loc[thiswav, key] = int(spd[key])
                        print(key, int(spd[key]) )
                    else:
                        df.loc[thiswav, key] = 0
                keymax = max(spd, key=spd.get)
                print('new_subclass = ',keymax)
                df.loc[thiswav, 'new_subclass']=keymax 
                checked = True
            if checked:
                df.loc[thiswav, 'checked']=True
        except: # BAD
            print('Input may have been faulty. Skipping event')
            #return empty_df, False
            return False
        else: # GOOD - RECLASSIFIED
            #df.drop(columns=['path'])                
            #return df, False     
            return False
    else: # BAD
        print("%s not found" % mseedpath)
        #return empty_df, False
        return False

#

def remove_marked_events(df):  
    dfsubset = df
    n_all = len(dfsubset.index)
    dfsubset = dfsubset[dfsubset['delete']==False]
    n_after_delete = len(dfsubset.index)
    dfsubset = dfsubset[dfsubset['ignore']==False]
    n_after_ignore = len(dfsubset.index)
    dfsubset = dfsubset[dfsubset['split']==False]
    n_after_split = len(dfsubset.index)
    n_delete = n_all - n_after_delete
    n_ignore = n_after_delete - n_after_ignore
    n_split = n_after_ignore - n_after_split
    print('Removed events: Marked to:')
    print('- split ', n_split)
    print('- delete ', n_delete)
    print('- ignore ', n_ignore)
    print('Catalog down from %d to %d events' % (n_all, n_after_split))
    print(' ')  
    return dfsubset

#

def _merge_dataframes(df_dict, accepted_subclasses):
    frames = []
    for subclass in accepted_subclasses:
        if subclass in df_dict.keys():
            frames.append(df_dict[subclass]) 
    return pd.concat(frames, sort=True) 

#

def _count_by_subclass(df):
    checked = df[df['checked']==True]
    unchecked = df[df['checked']==False]
    
    print('Event counts:')
    
    #print(checked.columns)
    
    if len(checked.index)>0:
        print('Checked events: %d' % len(checked.index) )
        checked_by_subclass = checked.groupby("new_subclass")
        print(checked_by_subclass['filetime'].count())
        
    if len(unchecked.index)>0:
        print('Unchecked events: %d' % len(unchecked.index) )
        unchecked_by_subclass = unchecked.groupby("new_subclass")
        print(unchecked_by_subclass['filetime'].count()) 

    print('Events by weight / quality threshold')
    print(df.groupby('weight').sfile.count())

#    
    
#def to_AAA(df, subclasses_for_ML, outfile, SEISAN_DATA, ignore_extra_columns=True, copy_to_path=None):
def to_AAA(df2, subclasses_for_ML, outfile, ignore_extra_columns=True):
    
    """
    create output file for AAA
    """ 
    df = df2.copy()
    included_subclasses = subclasses_for_ML.copy()
    #minweight = 0
    #df = df[df['weight']>=minweight]
    df = df[df['checked']==True]
    if len(df.index)<500:
        print('You need to reclassify some events first')
        return    
    print('Checked events: %d' % len(df.index)) 
    print(' ') 

    """
    print('Now we have the following number of events by subclass:')
    L = []
    for i, subclass in enumerate(subclasses_for_ML):
        df_subclass = df[df['new_subclass']==subclass]
        L.append(len(df_subclass.index))
        print('- %s: %d' % (subclass, L[i]))
    """
    
    print('all subclasses = ',df['new_subclass'].unique())
    print('included subclasses = ',included_subclasses)

    minnumevents = 10
    removed_subclasses = []
    for subclass in df['new_subclass'].unique():
        if subclass in included_subclasses:
            df_subclass = df[df['new_subclass']==subclass]
            if len(df_subclass.index) < minnumevents:
                print('Eliminating subclass %s' % subclass)
                included_subclasses.remove(subclass)
                removed_subclasses.append(subclass)
            
    print('The subclasses for machine learning are %s.' % ''.join(included_subclasses) )
    print('Removed subclasses %s' % ''.join(removed_subclasses) )
    df.sort_values(by='filetime',inplace=True)
    
    frames = []
    df_dict = {}
    for subclass in included_subclasses:
        df_dict[subclass] = df[df['new_subclass']==subclass]
    dfAAA = _merge_dataframes(df_dict, included_subclasses)    
    
    """
    df_list = []
    df.sort_values(by='filetime',inplace=True)
    for i, row in df.iterrows():
        row['new_subclass'] = row['new_subclass'].strip()
        if row['new_subclass'] in included_subclasses:
            df_list.append(row)
    df = pd.DataFrame(df_list)
    """
    
    print('Here is the FINAL list of events by subclass and whether they have been checked:')
    _count_by_subclass(dfAAA)    
    
    #df.rename(columns = {'twin':'duration'}, inplace = True)
    dfAAA.rename(columns = {'twin':'length'}, inplace = True)
    dfAAA['f0']=0.5
    dfAAA['f1']=25.0   
    
    """
    confidence_threshold = int(input('What is the minimum confidence percentage (e.g. 50) you wish to use ? This should be an integer.'))
    if confidence_threshold:
        good_indices = []
        for subclass in subclasses_for_ML:
            good_indices.expand(df[subclass] >= confidence_threshold)
        df = df.iloc[good_indices]
    print(df[['filetime', 'new_subclass']])
    """
    
    # subset and rename columns for output
    print(dfAAA.columns)
    if not 'path' in dfAAA.columns:
        dfAAA.reset_index(inplace=True)
    if ignore_extra_columns:
        dfAAA=dfAAA[['new_subclass','year','month','day','hour','minute','second','length','path']]     
    dfAAA.rename(columns = {'new_subclass':'class'}, inplace = True)
    
    """
    if copy_to_path:
        thisdir = os.getcwd()
        os.chdir(SEISAN_DATA)
        for i, row in dfAAA.iterrows():
            subclass = row['class']
            wavfile = row['path'] # relative path like ./WAV/DB/YYYY/MM/WAVfile
            newdir = os.path.join(copy_to_path, subclass)
            newfile = os.path.join(newdir, os.path.basename(wavfile))
            dfAAA.loc[i, 'path']=newfile
            if not os.path.exists(newdir):
                os.mkdir(newdir)
            if not os.path.exists(newfile):
                print('cp %s %s' % (wavfile, newfile) )
                os.system('cp %s %s' % (wavfile, newfile) )
        print('WAVfiles copied to ',copy_to_path)
        os.chdir(thisdir)
    """
    
    if os.path.exists(outfile):
        outfile_old = outfile + '.old'
        os.system("cp %s %s" % (outfile, outfile_old))
    _df2file_without_index(dfAAA, outfile)
    print('Catalog CSV saved to ',outfile)
    
#    
    
def report_checked_events(dfall, subclasses_for_ML):
    # Sanity check against AAA writer
    df = dfall.copy()
    df = df[df['checked']==True]
    if len(df.index)==0:
        print('You need to reclassify some events first')
        return
    print('total checked events = %d' % len(df.index))
    df = df[df['ignore']==False]
    df = df[df['delete']==False]
    df = df[df['split']==False]
    print('total classified events = %d' % len(df.index))
    frames = []
    for subclass in subclasses_for_ML:
        dfs = df[df['new_subclass']==subclass]
        print(subclass, len(dfs.index))
        frames.append(dfs)
    newdf = pd.concat(frames)

    L0 = len(newdf.index)
    print('total events matching ML subclasses = %d' % L0)

    differentdf = newdf[newdf['subclass']!=newdf['new_subclass']]
    L1 = len(differentdf.index)
    print('total reclassified events = %d' % L1)

    samedf = newdf[newdf['subclass']==newdf['new_subclass']]
    L2 = len(samedf.index)
    print('total already correctly classified events = %d' % L2)

    print('Error rate = %.1f%%' % (L1*100/L0))

    print(L0, L1+L2)

    #newdf.sort_values(by=['weight'],inplace=True,ascending=False)
    #print(newdf[['subclass', 'R', 'r', 'e', 'l', 'h', 't', 'new_subclass', 'weight']].to_string())

#
    
def _make_parent_dir(somefile):
    somedir = os.path.dirname(somefile)
    if not os.path.exists(somedir):
        os.makedirs(somedir)
    if not os.path.exists(somedir):
        print('Cannot create directory for %s' % somefile)
        exit()        

def _df2file_without_index(df, catfile, indexcol=None):
    df = df.reset_index()  
    if indexcol:
        if not indexcol in df.columns:
            df.rename(columns = {'index':indexcol})
    df.drop(df.filter(regex="Unname"),axis=1, inplace=True)
    if catfile[-3:]=='csv':
        df.to_csv(catfile, index=False)
    if catfile[-3:]=='pkl':
        df.to_pickle(catfile) 

def number_of_checked_events(df):
    df2 = df[df['checked']==True]
    print('Got %d checked event' % len(df2.index) )


# FUNCTIONS END ##########################################
##########################################################
# MAIN CODE BEGINS #######################################

# setup main directory paths & input/output files
SEISAN_DATA = os.path.join( os.getenv('HOME'),'DATA','MVO') # e.g. /home/user/seismo
#pandaSeisDir = os.path.join(SEISAN_DATA, 'pandaSeis') # e.g. /home/user/seismo/pandaSeis
pandaSeisDir = os.path.join(SEISAN_DATA, 'miniseed_c') # e.g. /home/user/seismo/pandaSeis
SEISAN_DB = 'MVOE_' # e.g. the seisan database name (e.g. MVOE_) under /home/user/seismo/WAV and /home/user/seismo/REA
pandaSeisDBDir = os.path.join(pandaSeisDir, SEISAN_DB) # e.g. /home/user/seismo/pandaSeis/MVOE_
AAA_DATA_DIR = os.path.join(SEISAN_DATA, 'MachineLearning', SEISAN_DB) # e.g. /home/user/seismo/MachineLearning/MVOE_
#master_event_catalog = os.path.join(AAA_DATA_DIR, 'labelling', '%s11_master_catalog.csv' % SEISAN_DB)
master_event_catalog = '/Users/thompsong/DATA/MVO/MachineLearning/MVOE_/labelling/11_merged_catalog.csv'
master_event_catalog_original = os.path.join(AAA_DATA_DIR, 'original', '%s11_master_catalog_original.csv' % SEISAN_DB)
aaa_input_file = os.path.join(AAA_DATA_DIR, 'runAAA', '%s11_labelled_events.csv' % SEISAN_DB) 
_make_parent_dir(master_event_catalog)
_make_parent_dir(master_event_catalog_original)
_make_parent_dir(aaa_input_file)
catalog_pickle_file = os.path.basename(master_event_catalog).replace('.csv','.pkl') 

# read items from Seisan software setup
volcanodefcsv = 'CSVfiles/volcano_def.csv'
if not os.path.exists(volcanodefcsv):
    print('%s not found. exiting' % volcanodefcsv)
    exit()
subclass_mapping = read_volcano_def(volcanodefcsv) # subclasses allowed for classification
seisan_subclasses = subclass_mapping['subclass'].values.tolist() # append('g') as needed, it is not an allowed subclass
#seisan_etypes = subclass_mapping['etype'].values.tolist()
subclasses_for_ML = ['D', 'R', 'r', 'e', 'l', 'h', 't'] # subclasses allowed for Machine Learning
station0hypfile = os.path.join(SEISAN_DATA, 'DAT', 'STATION0_MVO.HYP')

# Load a previously created master event catalog from AAA_DATA_DIR - or create a new one from pandaSeisDBDir
if os.path.exists(catalog_pickle_file):
    try:
        dfall = pd.read_pickle(catalog_pickle_file)   
    except:
        dfall = pd.read_csv(master_event_catalog) # how do i ignore the index?
elif os.path.exists(master_event_catalog): # load the one that exists from AAA_DATA_DIR, and trim the columns
    dfall = pd.read_csv(master_event_catalog) # how do i ignore the index?
else: # create one 
    print('%s does not exist, will try to create a new master event catalog from pandaSeisDBDir' % master_event_catalog)
    dfall = build_master_event_catalog(pandaSeisDBDir, SEISAN_DB, master_event_catalog_original, subclasses_for_ML)
    if os.path.exists(master_event_catalog_original):
        print(master_event_catalog_original, ' created')
        print('Copying to %s' % master_event_catalog)
        os.system("cp %s %s" % (master_event_catalog_original, master_event_catalog))
    else:
        print('Unable to make ',master_event_catalog_original, '. Quitting.')
        exit()
   
number_of_checked_events(dfall)
 
# ensure we always have the same index and columns
good_columns = []
for thiscol in dfall.columns:
    if 'ntitle' not in thiscol:
        good_columns.append(thiscol)
dfall = dfall[good_columns] # subset to correct columns
dfall.set_index('path', inplace=True) # try this
dfall.sort_index(inplace=True)   
    
    
number_of_checked_events(dfall)
# LOOPING BEGINS HERE ####################################

iterate_again = False # changed this back to do the loop
changes_made = True
while iterate_again:

    # get/update the fingerprints of each event class
    #fingerprints = get_weighted_fingerprints(dfall, subclasses_for_ML, N=100, exclude_checked=False)
    fingerprints = get_fingerprints(dfall, subclasses_for_ML, N=100, exclude_checked=False)
    save_fingerprints(fingerprints, subclasses_for_ML)
    
    # get index (path) of next event
    nextwav = select_next_event(dfall, subclasses_for_ML)
    print('Selected ',nextwav)
    
    # manually QC the next event. each time we choose the class with least checked examples
    bool_quit = qc_event(dfall, nextwav, subclasses_for_ML, seisan_subclasses, fingerprints, pandaSeisDBDir, station0hypfile)
    if bool_quit:
        iterate_again=False
    else: # changes made
        changes_made = True
        _df2file_without_index(dfall.copy(), catalog_pickle_file, 'path')     

# LOOPING ENDS HERE ######################################        


# Save
if changes_made:
    
    print('Quitting. Saving changes.')
    
    print('Updating %s' % master_event_catalog)
    _df2file_without_index(dfall, master_event_catalog, 'path')

    # remove events we marked for deletion, splitting or to ignore
    dfsubset = remove_marked_events(dfall)
    number_of_checked_events(dfsubset)

    #to_AAA(dfsubset, subclasses_for_ML, aaa_infile, SEISAN_DATA, ignore_extra_columns=False)
    print('Updating %s' % aaa_input_file)
    to_AAA(dfsubset, subclasses_for_ML, aaa_input_file, ignore_extra_columns=False) # need filetime column

    report_checked_events(dfall, subclasses_for_ML)
else:
    print('Quitting. No changes made.')
print('The end')


