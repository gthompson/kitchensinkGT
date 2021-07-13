#!/usr/bin/env python
#import obspy.core

import libMVO
#import trace_quality_control as qc
import os
import glob
import datetime as dt

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
            wavfile = mvocat.iloc[i].wavfilepath
            if wavfile.empty:
                continue
            try:
                os.path.exists(wavfile)
            except:
                print(wavfile[40:])
            try:
                if not os.path.exists(wavfile):
                    continue
            except:
                continue
            st = obspy.read(wavfile).select(id = traceid)
            if not len(st)==1:
                continue
            data = st[0].data
            if libMVO.check0andMinus1(data):
                quality_index[i] = 1
            else:
                pass
                #print("There are consecutive 0s or -1s in this recording")
        print('Done!')  
        mvocat['quality_index'] = quality_index
        mvocat.to_csv(mvocatqcfile)
        
    # add a fixed_trace_id column to dataframe - might want to make this first step
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
                libMVO.fix_trace_id(tr)
                fixed_ids[i]=tr.id
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

    # check trace ID for WAVfiles?
    days_df = pd.DataFrame()
    list_of_csv_files = sorted(glob.glob('./MVOE_wavfiles*.csv'))
    all_trace_ids = list()
    for csvfile in list_of_csv_files:
        df = pd.read_csv(csvfile)
        df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
    # loop over every day in this month & make a summary row for each channel
        old_trace_ids = df.traceid.unique()
        all_trace_ids = list(set().union(all_trace_ids, old_trace_ids))
    print('There are %d different traceids for the MVOE_ network' % len(all_trace_ids))
    for oldid in all_trace_ids:
        if oldid.find('MB')>-1:
            newid = libMVO.correct_nslc(oldid, 75.19)  
        else:
            newid = libMVO.correct_nslc(oldid, 100.0)
        print(oldid, '->', newid)             
                
            
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
          

    """ From count_traces_per_day_MVOE.py """
    days_df = pd.DataFrame()
    # Load each monthly catalog
    #list_of_csv_files = sorted(glob.glob('./ASNE_wavfiles*.csv'))
    list_of_csv_files = sorted(glob.glob('./MVOE_wavfiles*.csv'))
    print(list_of_csv_files)
    all_trace_ids = list()
    # get a list of all traceids
    for csvfile in list_of_csv_files:
        print(csvfile)
        df = pd.read_csv(csvfile)
        #print(df)
        df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
        # loop over every day in this month & make a summary row for each channel
        new_trace_ids = df.fixedid.unique()
        all_trace_ids = list(set().union(all_trace_ids, new_trace_ids))   
    print(all_trace_ids)

    dailytraceid_df = pd.DataFrame(columns=['yyyymmdd', 'traceid', 'count']) 
    for csvfile in list_of_csv_files:
        print(csvfile)
        df = pd.read_csv(csvfile)
        df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
        yyyy = int(csvfile[-10:-6])
        mm = int(csvfile[-6:-4])
        start_date = dt.date(yyyy, mm, 1)
        emm = mm+1
        eyyy=yyyy
        if emm==13:
            emm=1
            eyyy = yyyy+1
        end_date = dt.date(eyyy, emm, 1)
        delta = dt.timedelta(days=1)
        while start_date < end_date: # loop over all days in this month
            yyyymmdd = start_date.strftime("%Y%m%d")
            print(yyyymmdd)
            sfiledatetime = pd.to_datetime(df.datetime)
            daycat = df.loc[sfiledatetime.dt.date == start_date] # get rows just for this day
            #daytraceid_hash = daycat.traceid.value_counts() 
            print(daycat)
            print(daycat.shape)
            for thistraceid in all_trace_ids:
                if not daycat.empty: 
                    print(thistraceid)
                    try:
                        thisdaytrace_df = daycat.loc[daycat['fixedid'] == thistraceid]
                        thiscount = thisdaytrace_df.shape[0] 
                        #print('*** it works ***')
                    except:
                        print('thistraceid = %s' % thistraceid)
                        print(daycat)
                        barf 
                else:
                    thiscount = 0
                #thisrow = pd.DataFrame(yyyymmdd, thistraceid, thiscount)
                thisrow = {'yyyymmdd':yyyymmdd, 'traceid':thistraceid, 'count':thiscount}
                dailytraceid_df = dailytraceid_df.append(thisrow, ignore_index=True)
                #print(thisrow)
                
            start_date += delta # add a day  
        #print('Done with this csvfile')
    dailytraceid_df.to_csv('MVOE_dailytraceid_wavfiles_df.csv')  

            
    """ See station_ontimes_MVOE.ipynb. This might already be in the ASN paper repo """

    
    """ from csv2qc.py """
    list_of_csv_files = sorted(glob.glob('ASNE_catalog??????.csv'))
    #frames = []
    for csvfile in list_of_csv_files:
        df = pd.read_csv(csvfile)
        #frames.append(df)

        # fix column names
        df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')

        nrows, ncolumns = df.shape
        print(nrows)
        quality_index = np.zeros(nrows)
        df['quality_index'] = quality_index
        fixed_ids = df.traceid.tolist() # a list of all trace ids in the csvfile

        for i in range(nrows):
            if i % 100 == 0:
                print('Done %d of %d' % ( i, nrows, ))

            # code from fix_montserrat_traceids.ipynb
            traceid = fixed_ids[i]
            wavfilepath = mvocat.iloc[i].wavfilepath
            #print(wavfilepath)
            if os.path.exists(wavfilepath):
                st = obspy.read(wavfilepath)
                tr = st[df.iloc[i].tracenum]
                print(tr)
                tr = qc.compute_metrics(tr)
                tr2 = fix.fix_trace_ids(tr)
                fixed_ids[i]=tr2.id # fix the list
            else:
                print('does not exist')

            # code from qc_traces_from_CSV.ipynb
            #traceid = mvocat.iloc[i].traceid
            #wavfile = mvocat.iloc[i].wavfilepath
            #if wavfile.empty:
            #    continue
            #try:
            #    os.path.exists(wavfile)
            #except:
            #    print(wavfile[40:])
            #try:
            #    if not os.path.exists(wavfile):
            #        continue
            #except:
            #    continue
            #st = obspy.read(wavfile).select(id = traceid)
            if not len(st)==1:
                continue
            data = st[0].data
            if check0andMinus1(data):
                quality_index[i] = 1
            else:
                pass
                #print("There is consecutive 0 or -1. in this recording")
        print('Done!')
        df['quality_index'] = quality_index
        qc_csvfile = csvfile[0:12] + "_qc_.csv"
        df.to_csv(qc_csvfile)

