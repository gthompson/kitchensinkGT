#!/usr/bin/env python
#import obspy.core

"""
## Step 1: setup a .bashrc_montserrat file in your ~ (HOME) directory
### Example of file contents:
export PYTHONPATH="/Users/thompsong/src/Alexis_Montserrat_codes"
export PATH=$PATH:/Users/thompsong/src/Alexis_Montserrat_codes
export SEISAN_DATA="/Users/thompsong/Desktop/seismo"

## Step 2: source this from your .bashrc by doing
echo "source ~/.bashrc_montserrat" >> ~/.bashrc

## Step 3: start a new bash terminal

## Step 4: run seisandb2csv.py
cd $SEISAN_DATA
seisandb2csv.py

## Step 5: A CSV file like MVOE_catalog200501.csv should now have been created for each month that S-files are in your Seisan database, e.g. under $SEISAN_DATA/REA. You could now use standard Linux commands, e.g.

## (a) grep out 1 channel only
grep MBWH.Z.BH MVOE_catalog*.csv > MBWH.BHZ.csv

## (b) grep out 1 channel and count duration 
grep MBWH.Z.BH MVOE_catalog*.csv | awk -F',' '{print $4}' 

## (c) count the number of each volcanic subclass
awk -F',' '{print $3}' MVOE_*.csv | sort | uniq -c

## (d) count the number of seconds total for each different channel
awk -F',' '{a[$9]+=$4}END{for(i in a) print i,a[i]}' MV*csv | sort

## Step 6: CSV files can directly be imported in Python dataframes and manipulated there.

## Step 7: Alexis, this would be a great place to also put any *.pynb and *.json and other *.py codes related to running malfante/AAA codes on the Montserrat data

Glenn Thompson 2019/10/16, at IPGP

"""

def check0andMinus1(liste):
    liste=list(liste)
    listStr=''.join(str(i) for i in liste)
    if  "000000000000" in listStr or "-1-1-1-1-1-1-1-1" in listStr :
        return False
    else:
        return True

def fix_trace_ids(tr): # add stuff here to correct sampling rate too
    # fix the network, channel and location
    network = 'MV'
    tr.stats['network']=network
    sta = tr.stats['station'].strip()
    chan = tr.stats['channel'].strip()
    if chan=='PRS' or chan=='APS':
        chan='BDF'
    else:
        if chan[0]=='A':
            if tr.stats['location'] == 'J':
                bandcode = 'S'
            #else:
                #bandcode = 'B'
        else:
            if chan[1]=='B':
                bandcode = 'B'
            else:
                bandcode = chan[0]
            instrumentcode = 'H'
            if len(chan)==2:
                orientationcode = tr.stats['location']
            else:
                orientationcode = chan[2]
            chan = bandcode + instrumentcode + orientationcode

    if chan[0]=='A':
        print(tr.stats)
        print(chan)
        sys.exit('bugger!')
    tr.stats['channel'] = chan
    tr.stats['location']='--'
    return tr

if __name__ == '__main__':
    import os, sys, glob
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    cwd = os.getcwd()  
    SEISAN_DB = 'MVOE_'
    SEISAN_DATA = os.environ['SEISAN_DATA']
    CSVpath = os.path.dirname(SEISAN_DATA)
    os.chdir(CSVpath)
    
    """ #####################################################
                      PROCESS REA DIRECTORY 
        ##################################################### """
    
    """ 
    Traverse a Seisan database structure by REA/year/month, find all the S-files and corresponding WAV files, 
    and generate a summary line for each trace. 
    
    Will generate nothing for an Sfile for which no corresponding WAVfile found 
    
    Output files are like DBNAMEcatalogYYYYMM.csv
    
    This is a wrapper to run the "generate_monthly_csv" method of the Seisan_Catalog class
    
    The CSV files are written into the parent directory of the path stored in the SEISAN_DATA environment variable. 
    SEISAN_DATA is the parent of the WAV and REA directories, e.g. /Users/USERNAME/montserrat_data/seismo 
    In this case, the CSV files would be written to /Users/USERNAME/montserrat_data/ 
    
    The default database is "MVOE_", so the CSV files are named like MVOE_catalogYYYYMM.CSV
    
    If a CSV already exists, it will not be overwritten. 
    So if you want them to be overwritten, delete them before running this
    """
    monthly_csv_files = sorted(glob.glob('MVOE_catalog*.csv'))
    if len(monthly_csv_files)==0:
        yyyylist = sorted(glob.glob(os.path.join(SEISAN_DATA, 'REA', SEISAN_DB, '????')))
        print(yyyylist)
        flag_sfiles_only = True
        #flag_sfiles_only = False
        for yyyy in yyyylist:
            mmlist = sorted(glob.glob(os.path.join(yyyy, '??')))
            for mm in mmlist:
                Seisan_Catalog.generate_monthly_csv(mm, flag_sfiles_only)    
    monthly_csv_files = sorted(glob.glob('MVOE_catalog*.csv'))
    
    """ read the monthly catalog files into a single dataframe """
    mvocatfile = 'MVOE_CATALOG.csv'
    if os.path.exist(mvocatfile):
        mvocat = pd.read_csv(mvocatfile)
    else:
        frames = []
        for monthcsvfile in monthly_csv_files:
            thisframe = pd.read_csv(monthcsvfile)
            frames.append(thisframe)
        mvocat=pd.concat(frames)
        mvocat.reset_index()
        #mvocat.columns = mvocat.columns.str.lstrip()
        mvocat.columns = mvocat.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
        mvocat.to_csv(mvocatfile)
  
    """ add a quality_index column to dataframe """
    mvocatqcfile = 'MVOE_CATALOG_QC.csv'
    if os.path.exist(mvocatqcfile):
        mvocat = pd.read_csv(mvoqccatfile)
    else:    
        nrows, ncolumns = mvocat.shape 
        quality_index = np.zeros(nrows)
        mvocat['quality_index'] = quality_index
        for i in range(nrows):
            print(i, end = ' ')
            traceid = mvocat.iloc[i].traceid
            st = obspy.read(mvocat.iloc[i].wavfilepath).select(id = traceid)
            data = st[0].data
            if check0andMinus1(data):
                quality_index[i] = 1
            else:
                pass
                #print("There are consecutive 0s or -1s in this recording")
        print('Done!')  
        mvocat['quality_index'] = quality_index
        mvocat.to_csv(mvocatqcfile)
  
    # add a fixed_trace_id column to dataframe
    mvocatfixedIDfile = 'MVOE_CATALOG_QC_fixedID.csv'
    if os.path.exist(mvocatfixedIDfile):
        mvocat = pd.read_csv(mvocatfixedIDfile)
    else:
        nrows, ncolumns = mvocat.shape
        mvocat.columns = mvocat.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
        fixed_ids = mvocat.traceid.tolist()
        for i in range(nrows):
            if i % 100 == 0:
                print('Done %d of %d' % ( i, nrows, ))
            traceid = fixed_ids[i]
            wavfilepath = mvocat.iloc[i].wavfilepath
            print(wavfilepath)
            if os.path.exists(wavfilepath):
                st = obspy.read(wavfilepath)
                tr = st[mvocat.iloc[i].tracenum]
                print(tr)
                tr2 = fix.fix_trace_ids(tr)
                fixed_ids[i]=tr2.id
            else:
                print('does not exist')
        mvocat['fixed_trace_id'] = fixed_ids
        mvocat.to_csv(mvocatfixedIDfile)
    
    # various counts
    mvocat.mainclass.value_counts()
    mvocat.subclass.value_counts()    
    LVcat = mvocat.loc[mvocat.mainclass == 'LV']
    LVcat.subclass.value_counts()
    len(LVcat.wavfilepath.unique())     
    
    """ #####################################################
                      PROCESS WAV DIRECTORY 
        ##################################################### """  
    
    """ 
    Traverse a Seisan database structure by WAV/db/year/month finding all the WAV-files
    Generate a summary line for each trace
    Output files are like DBNAMEwavfilesYYYYMM.csv
    These will later be read into Python dataframes
    """
    yyyylist = sorted(glob.glob(os.path.join(SEISAN_DATA, 'WAV', SEISAN_DB, '????')))
    print(yyyylist)
    for yyyy in yyyylist:
        mmlist = sorted(glob.glob(os.path.join(yyyy, '??')))
        for mm in mmlist:
            #print("Processing %s:" % mm)
            Seisan_Catalog.generate_monthly_wav_csv(mm)    

            
            
    # data completeness by station
    oncsvfile = 'station_onness_SIMPLE.csv' # have not found code that generated this
    if os.path.exists(oncsvfile):
        df = pd.read_csv(oncsvfile)
        df = df.fillna(0)
        df = df.replace('1',1)
        data = df.to_numpy()
        data = np.nan_to_num(data)
        print(data)
        data.shape
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(data, aspect='auto', cmap=plt.cm.gray, interpolation='nearest')
        xticklabels = ('MBBE','MBBY','MBFL','MBFR','MBGA','MBGB','MBGE','MBGH','MBHA','MBLG','MBLY','MBMH','MBRV','MBRY','MBSS','MBWH')
        plt.xticks(np.arange(16), xticklabels)
        ax.set_xticklabels(xticklabels, rotation = 90)
        yticklabels = ('1996/10','1997/01','1998/01','1999/01','2000/01','2001/01','2002/01','2003/01','2004/01','2005/01','2006/01','2007/01')
        plt.yticks([0, 3, 15, 27, 39, 51, 63, 75, 87, 99, 111, 123, 135], yticklabels)
        plt.tight_layout()
        plt.savefig('station_onness.png',dpi=200)
          
