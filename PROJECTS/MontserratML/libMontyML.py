import os
from glob import glob
import numpy as np
import pandas as pd
from IPython.display import clear_output
from statsmodels.stats.weightstats import DescrStatsW
from obspy.core import read
import matplotlib.pyplot as plt

def _select_next_event(dfall, subclasses_for_ML):
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
        print('next event is of subclass = %s' % nextclass)

        is_subclass = unchecked['subclass']==nextclass
        df = unchecked[is_subclass]

        # mechanism to weight out checked events
        df = df[df['ignore']==False]

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
            df = df.head(1)
            if 'bandratio_[1.0_6.0_11.0]' in df.columns:
                df.rename(columns = {'bandratio_[1.0_6.0_11.0]':'band_ratio'}, inplace = True)
            return df
        else:
            subclasses.remove(nextclass)
    return None

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



def _guess_subclass(row, fingerprints, subclasses_for_ML):
    chance_of = {}
    for item in subclasses_for_ML:
        chance_of[item]=0.0
    
    params = ['peaktime', 'kurtosis', 'medianF', 'peakF', 'bw_min', 'bw_max', 'band_ratio']
    for subclass in fingerprints.keys():
        fingerprint = fingerprints[subclass]
        #fingerprint.reset_index(inplace=True)
        #print(fingerprint.columns)
        for param in params:
            thisval = row[param]
            
            # test against mean+/-std
            meanval = fingerprint[param]['mean']
            stdval = fingerprint[param]['std']
            minus1sigma = meanval - stdval
            plus1sigma = meanval + stdval
            
            if thisval > minus1sigma and thisval < plus1sigma:
                weight = 1.0 - abs(thisval-meanval)/stdval
                chance_of[subclass] += weight
                
            # test against 25-75% percentile
            medianval = fingerprint[param]['50%']
            val25 = fingerprint[param]['25%']
            val75 = fingerprint[param]['75%']
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


def qc_event(dfall, subclasses_for_ML, seisan_subclasses, fingerprints, SEISAN_DATA, station_locationsDF):
    df = _select_next_event(dfall, subclasses_for_ML)
    subclass = df.iloc[0]['new_subclass']

    for index, row in df.iterrows():
        clear_output(wait=True)
        picklepath = os.path.join(SEISAN_DATA, row.path.replace('./WAV','PICKLE') + '.pickle')
        
        if os.path.exists(picklepath):
            
            # plot whole event file
            st = read(picklepath) #.select(station='MBWH')
            st.filter('bandpass', freqmin=0.5, freqmax=25.0, corners=4)
            deconvolve_instrument_response(st)

            # compute ampeng to show later
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
            
            # plot all
            st.plot(equal_scale=False);

            # also plot a fixed 30-seconds around the peaktime
            starttime = st[0].stats.starttime + row['peaktime']-10
            endtime = starttime + 20
            st.trim(starttime=starttime, endtime=endtime, pad=True, fill_value=None)
            st.plot(equal_scale=False)                  

            # Show a map of amplitude distribution
            try:
                add_station_locations(st, station_locationsDF)
                plot_amplitude_locations(st)
            except:
                pass

            # TEXT OUTPUT STARTS HERE
            print(' ')
            print('Checked events currently:')
            #print(dfall[dfall['checked']==True].groupby('new_subclass').size())
            print(row)
            print(' ')
            print('Loaded %s' % picklepath)
            
            # station amplitudes and energies
            #print(dfenergy)

            # trace df metrics
            csvpath = picklepath.replace('.pickle', '.csv')
            tracedf = pd.read_csv(csvpath)
            if 'bandratio_[1.0_6.0_11.0]' in tracedf.columns:
                tracedf.rename(columns = {'bandratio_[1.0_6.0_11.0]':'band_ratio'}, inplace = True)
            tracedf.sort_values(by=['peakamp'], ascending=False, inplace=True)
            print(' ')
            print(tracedf[['id', 'medianF', 'bw_min', 'peakF', 'bw_max', 'band_ratio', 'kurtosis']])

            # fingerprint guesses
            print(' ')
            _guess_subclass(row, fingerprints, subclasses_for_ML)
            print(suggested_weight)

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
                if new_subclass == 'q': 
                    return df, True
                if new_subclass == 's':
                      df.loc[index, 'split'] = True
                      checked = True
                if new_subclass == 'i':
                      df.loc[index, 'ignore'] = True 
                      checked = True
                if new_subclass == 'd':
                      df.loc[index, 'delete'] = True 
                      checked = True
                if not checked:                         
                    if not ',' in new_subclass: # convert to a subclass, percentage string
                        new_subclass = new_subclass + ', 100'
                    spl = new_subclass.split(',') # split string to subclass probability list 
                    if len(spl) % 2 == 1:
                        df.loc[index, 'weight'] = int(spl.pop())
                    spd = {spl[a].strip():spl[a + 1] for a in range(0, len(spl), 2)} # subclass probability dict
                    for key in subclasses_for_ML:
                        if key in spd.keys():
                            df.loc[index, key] = int(spd[key])
                            print(key, int(spd[key]) )
                        else:
                            df.loc[index, key] = 0
                    keymax = max(spd, key=spd.get)
                    print('new_subclass = ',keymax)
                    df.loc[index, 'new_subclass']=keymax  
                    checked = True
                if checked:
                    df.loc[index, 'checked']=True
            except:
                print('Input may have been faulty. Skipping event')
                return False
                    
    return df, False




def parse_STATION0HYP(station0hypfile):
    """ file sample
    012345678901234567890123456 
      MWNH1644.53N 6211.45W 407
      MWNL1644.53N 6211.45W 407
      MWZH1644.53N 6211.45W 407
      MWZL1644.53N 6211.45W 407
      0123456789012345678901234 
    """
    station_locations = []
    if os.path.exists(station0hypfile):
        fptr = open(station0hypfile,'r')
        for line in fptr:
            line = line.strip()
            if len(line)==25:
                station = line[0:4]
                latdeg = line[4:6]
                if latdeg != '16':
                    print(line)
                    continue
                latmin = line[6:8]
                latsec = line[9:11]
                hemisphere = line[11]
                if hemisphere == 'N':
                    latsign = 1
                elif hemisphere == 'S':
                    latsign = -1 
                else:
                    print(line)
                    continue
                londeg = line[13:15]
                lonmin = line[15:17]
                lonsec = line[18:20]
                lonsign = 1
                if line[20]=='W':
                    lonsign = -1
                elev = line[21:25].strip()
                station_dict = {}
                station_dict['name'] = station
                station_dict['lat'] = (float(latdeg) + float(latmin)/60 + float(latsec)/3600) * latsign
                station_dict['lon'] = (float(londeg) + float(lonmin)/60 + float(lonsec)/3600) * lonsign
                station_dict['elev'] = float(elev)
                
                station_locations.append(station_dict)
        fptr.close()
        return pd.DataFrame(station_locations)

def add_station_locations(st, station_locationsDF):
    for tr in st:
        df = station_locationsDF[station_locationsDF['name'] == tr.stats.station]
        df = df.reset_index()
        if len(df.index)==1:
            row = df.iloc[0]
            tr.stats['lat'] = row['lat']
            tr.stats['lon'] = row['lon']
            tr.stats['elev'] = row['elev']
        else:
            print('Found %d matching stations' % len(index))
            
def plot_amplitude_locations(st):
    plt.figure()
    x = []
    y = []
    s = []
    for tr in st:
        if 'lon' in tr.stats:
            x.append(tr.stats.lon)
            y.append(tr.stats.lat)
            s.append(np.max(np.absolute(tr.data)))
    if len(s)>0:
        maxs = np.max(s)
        s = s/maxs * 50
        plt.scatter(x, y, s, facecolors='none', edgecolors='r');
        for i, label in enumerate(tr.stats.station):
            #if tr.stats.channel[-1]=='Z':
            plt.text(x[i]*1.01, y[i]*1.01, label)

        volcano_dict = {}
        volcano_dict['name']='VOLC'
        volcano_dict['lat']=16.7166638
        volcano_dict['lon']=-62.1833326
        volcano_dict['elev']=1000
        plt.scatter(volcano_dict['lon'], volcano_dict['lat'], 300, marker='*')
        plt.show();

def deconvolve_instrument_response(st):
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
                this_inv = load_mvo_inventory(tr, '/Users/thompsong/DATA/MVO/CAL')
        if this_inv and tr.stats['units'] == 'Counts':
            tr.remove_response(inventory=this_inv, output="VEL")
            tr.stats['units'] = 'm/s'
            add_to_trace_history(tr, 'deconvolved')
        #else:
        #    st.remove(tr)
    

def read_volcano_def():
    filepath = './volcano_def.csv'
    subclass_df = pd.read_csv(filepath)
    subclass_df.columns = subclass_df.columns.str.lstrip()
    return subclass_df

def build_master_trace_catalog(SEISAN_PATH, catalogdf, subclass=None):
    frames = []
    cat = catalogdf.copy()
    if subclass:
        cat = cat[cat.subclass == subclass]
    for i,row in cat.iterrows():
        traceCSV = row['path'].replace('./WAV', SEISAN_PATH + '/PICKLE').replace('.pickle','') + '.csv'
        df = pd.read_csv(traceCSV)
        frames.append(df)
    traceDF = pd.concat(frames, sort=True)
    return traceDF
    
def build_master_event_catalog(csvdir, seisandbname, catalogfile, subclasses_for_ML, max_duration = 300):
    # load all the year/month CSV files
    csvfiles = glob(os.path.join(csvdir, 'reawav_%s??????.csv' % seisandbname))
    frames = []
    for csvfile in csvfiles:
        df = pd.read_csv(csvfile)
        frames.append(df) 
    dfall = pd.concat(frames, sort=True)
    dfall.set_index('filetime', inplace=True) # we will need this later to remerge
    dfall.sort_index(inplace=True)
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
    dfall.to_csv(catalogfile)
    
    return dfall

def _count_by_subclass(df):
    checked = df[df['checked']==True]
    unchecked = df[df['checked']==False]
    
    print('Event counts:')
    
    if len(checked.index)>0:
        print('Checked events: %d' % len(checked.index) )
        checked_by_subclass = checked.groupby("new_subclass")
        print(checked_by_subclass['path'].count())
        
    if len(unchecked.index)>0:
        print('Unchecked events: %d' % len(unchecked.index) )
        unchecked_by_subclass = unchecked.groupby("new_subclass")
        print(unchecked_by_subclass['path'].count()) 

    print('Events by weight / quality threshold')
    print(df.groupby('weight')['path'].count())
    
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
                if 'bandratio_[1.0_6.0_11.0]' in df_subclass.columns:
                    df_subclass.rename(columns = {'bandratio_[1.0_6.0_11.0]':'band_ratio'}, inplace = True)
                best_events_dict[subclass]=df_subclass

    return best_events_dict


def get_fingerprints(dfall, allowed_subclasses, N=100, exclude_checked=False):
    """
    All we do right now is a dataframe describe, so we return the stats of each column.
    
    wdf = DescrStatsW(df.x, weights=df.wt, ddof=1) 
    
    """
    #df = best_events_dict.groupby("subclass")
    fingerprints = {}
    best_events_dict = _select_best_events(dfall, allowed_subclasses, N=N, exclude_checked=exclude_checked)
    for subclass in allowed_subclasses:
        if subclass in best_events_dict.keys(): 
            print('Computing fingerprint for subclass ',subclass)
            #df_subclass = df.get_group(subclass)
            df_subclass = best_events_dict[subclass]     
            fingerprints[subclass] = df_subclass[[ 'peaktime', 
                'kurtosis', 'medianF', 'peakF', 'bw_min', 'bw_max', 'band_ratio']].describe()
            #print(fingerprints[subclass])
        else:
            fingerprints[subclass]=pd.read_csv('fingerprint_%s.csv' % subclass)
    return fingerprints

def save_fingerprints(fingerprints, allowed_subclasses):
    for subclass in allowed_subclasses:
        if subclass in fingerprints.keys(): 
            fingerprints[subclass].to_csv('fingerprint_%s.csv' % subclass)
      
def _merge_dataframes(df_dict, accepted_subclasses):
    frames = []
    for subclass in accepted_subclasses:
        if subclass in df_dict.keys():
            frames.append(df_dict[subclass]) 
    return pd.concat(frames, sort=True) 

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

def to_AAA(df, subclasses_for_ML, outfile, SEISAN_DATA, ignore_extra_columns=True, copy_to_path=None):
    """
    create output file for AAA
    """ 
    included_subclasses = subclasses_for_ML.copy()
    minweight = 0
    df = df[df['weight']>=minweight]
    df = df[df['checked']==True]
    
    print(' ') 
    print('Now we have the following number of events by subclass:')
    L = []
    for i, subclass in enumerate(subclasses_for_ML):
        df_subclass = df[df['new_subclass']==subclass]
        L.append(len(df_subclass.index))
        print('- %s: %d' % (subclass, L[i]))
    
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
    if ignore_extra_columns:
        dfAAA=dfAAA[['new_subclass','year','month','day','hour','minute','second','length','path']]     
    dfAAA.rename(columns = {'new_subclass':'class'}, inplace = True)
    
    
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
            
    dfAAA.to_csv(outfile)
    print('Catalog CSV saved to ',outfile)
    
def report_checked_events(dfall, subclasses_for_ML):
    # Sanity check against AAA writer
    df = dfall.copy()
    df = df[df['checked']==True]
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



