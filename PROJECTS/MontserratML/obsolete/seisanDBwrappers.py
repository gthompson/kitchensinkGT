#!/usr/bin/env python
import os, sys
from glob import glob
import numpy as np
import pandas as pd
import datetime as dt
from pprint import pprint
import matplotlib.pyplot as plt

from obspy import read, read_inventory, Stream
from obspy.io.xseed.core import _read_resp
from obspy.imaging.cm import obspy_sequential

LIBpath = os.path.join( os.getenv('HOME'),'src','kitchensinkGT', 'LIB')
sys.path.append(LIBpath)
from libMVO import fix_trace_id, inventory_fix_id_mvo
from metrics import process_trace, choose_best_traces, select_by_index_list, ampengfft
from libseisGT import Stream_min_starttime, detect_network_event
from Sfile import spath2datetime, Sfile #, printEvents

sys.path.append(os.path.join( os.getenv('HOME'),'src', 'icewebPy') )
import IceWeb

""" 
A set of tools to process Seisan databases (REA and WAV directories)
Created for Montserrat catalog work.
Note there are separate classes for loading individual Sfiles, AEFfiles and WAVfiles.
"""

def get_sfile_list(SEISAN_DATA, DB, startdate, enddate): 
    """
    make a list of Sfiles between 2 dates
    """

    event_list=[]
    reapath = os.path.join(SEISAN_DATA, 'REA', DB)
    years=list(range(startdate.year,enddate.year+1))
    for year in years:
        if year==enddate.year and year==startdate.year:
            months=list(range(startdate.month,enddate.month+1))
        elif year==startdate.year:
            months=list(range(startdate.month,13))
        elif year==enddate.year:
            months=list(range(1,enddate.month+1))
        else:
            months=list(range(1,13))
        for month in months:
            #print month
            yearmonthdir=os.path.join(reapath, "%04d" % year, "%02d" % month)
            flist=sorted(glob(os.path.join(yearmonthdir,"*L.S*")))
            for f in flist:
                #fdt = sfilename2datetime(f)
                fdt = spath2datetime(f)
                #print(f, fdt)
                if fdt>=startdate and fdt<enddate:
                    event_list.append(f)
    return event_list 


def generate_monthly_csv(mm, flag_sfiles_only=False):
# function generate_monthly_csv(monthdirs, flag_sfiles_only)
    print("Reading %s" % mm)
    sfileslist = sorted(glob(os.path.join(mm, "*")))
    if flag_sfiles_only:
        outfile = mm[-13:-8] + 'sfiles' + mm[-7:-3] + mm[-2:] + '.csv'
        hypofile = mm[-13:-8] + 'hypo' + mm[-7:-3] + mm[-2:] + '.csv'
    else:
        outfile = mm[-13:-8] + 'catalog' + mm[-7:-3] + mm[-2:] + '.csv'
    print("...Creating %s" % outfile)
    if os.path.exists(outfile):
        print("CSV file %s already exists. Not overwriting. If you want to overwrite then please delete from the command line" % outfile)
        return
    fptr = open(outfile,'w')
    if flag_sfiles_only:
        fptr.write('datetime,mainclass,subclass,sfilepath,analyst,wavfile1,wavfile2,numarrivals,nummagnitudes,numaefrows,duration\n') 
    else:
        fptr.write('datetime,mainclass,subclass,duration,wavfilepath,sampling_rate,npts,traceNum,traceID,sfilepath,analyst,fixedID,quality_factor,snr,highval,lowval\n') 
    for thissfile in sfileslist:
        s = Sfile(thissfile)
        if not np.isnan(s.latitude) and flag_sfiles_only:
            if not s.longitude:
                s.longitude = float('nan')
            if not s.depth:
                s.depth = float('nan')
            # write to hypo file
            hfptr = open(hypofile,'a+')
            if not os.path.exists(hypofile):
                hfptr.write('datetime,mainclass,subclass,latitude,longitude,depth,magnitude,magnitude_type,sfilepath\n') 
            if len(s.magnitude)>0:
                magnitude = s.magnitude[0]
                magnitude_type = s.magnitude_type[0]
            else:
                magnitude = float('nan')
                magnitude_type = ""
            try:
                hfptr.write("%s,%s,%s,%9.4f,%9.4f,%6.1f,%4.1f,%s,%s\n" % (s.filetime, s.mainclass, \
                    s.subclass, \
                    s.latitude, s.longitude, s.depth, magnitude, magnitude_type, \
                    thissfile ))
            except:
                print(s.latitude)
                print(s.longitude)
                print(s.depth)
                print(magnitude)
                barf
            hfptr.close()
        numwavfiles = len(s.wavfiles)
        tmp = thissfile.split('REA')
        shortsfile = 'REA' + tmp[1] 
        if not s.subclass:
            s.subclass = "_"
        if flag_sfiles_only:
            wavfile1 = "-"
            wavfile2 = "-"
            if numwavfiles>0:
                tmp = s.wavfiles[0].path.split('WAV')
                wavfile1 = 'WAV' + tmp[1]
            if numwavfiles>1:
                tmp = s.wavfiles[1].path.split('WAV')
                wavfile2 = 'WAV' + tmp[1]
            if len(s.aeffiles)>0:
                duration = s.aeffiles[0].trigger_window
            else:
                duration = 0.0
            if not duration:
                duration = 0.0

            try:
                fptr.write("%s,%s,%s,%s,%s,%s,%s,%2d,%1d,%2d,%5.1f\n" % (s.filetime, s.mainclass, \
                    s.subclass, \
                    thissfile, s.analyst, wavfile1, wavfile2, \
                    len(s.arrivals), len(s.magnitude), len(s.aefrows), duration  ))
            except:
                print(len(s.arrivals))
                print(len(s.magnitude))
                print(len(s.aefrows))
                print(duration)
                barf

        else:
            for thiswavfile in s.wavfiles:
                if not os.path.exists(thiswavfile.path):
                    continue
                try:
                    st = read(thiswavfile.path)
                except:
                    print("Processing %s: Cannot read wavfile %s" % (thissfile, thiswavfile.path,) )
                    continue
                tracenum = 0
                tmp = thiswavfile.path.split('WAV')
                wavfile1 = 'WAV' + tmp[1]
                for tr in st:
                     tr2, quality_factor, snr = qc.compute_metrics(tr)
                     duration = tr.stats.npts / tr.stats.sampling_rate
                     fptr.write("%s,%s,%s,%8.3f,%s,%3d,%6d,%02d,%s,%s,%s,%s,%d,%.2f,%.2f,%.2f\n" % (s.filetime, s.mainclass, \
                         s.subclass, duration, \
                         wavfile1, \
                         tr.stats.sampling_rate, tr.stats.npts, \
                         tracenum, tr.id, shortsfile, s.analyst, tr2.id, quality_factor, snr[0], snr[1], snr[2] ))
                     tracenum += 1
            
    fptr.close()


def sfilecsv_daycount(list_of_csv_files):
# function sfilecsv_daycount(list_of_csv_files)

    # Combine all the CSV files, and then summarize
    print('\n*************************')
    print('**** Overall summary ****')
    frames = []
    for csvfile in list_of_csv_files:
        df = pd.read_csv(csvfile)
        df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
        if not isinstance(frames, pd.DataFrame):
            frames = df
        else:
            df = df.reset_index(drop=True)
        frames = pd.concat([frames, df], axis=0)

    # how many of each mainclass?
    class_hash = frames.mainclass.value_counts()
    for ck in class_hash.keys():
        print(ck,class_hash[ck], end=' ')
    print(' ')

    # how many of each subclass?
    subclass_hash = frames.subclass.value_counts()
    for sk in subclass_hash.keys():
        print(sk,subclass_hash[sk], end=' ')
    print(' ')

    # how many by each analyst?
    analyst_hash = frames.analyst.value_counts()
    for ak in analyst_hash.keys():
        print(ak,analyst_hash[ak], end=' ')
    print(' ')

    print('***********************\n')

    # Now create the day summary dataframes
    # What we really want to do here is build a CSV file that contains date, and number of each mainclass and subclass for that date
    # use mainclass and subclass from above
    
    # create a new dataframe to do daily counts of eventtype
    cols = ['date']
    for ck in class_hash.keys():
        cols.append(ck)
    for sk in subclass_hash.keys():
        cols.append(sk)
    #eventtype_df = pd.DataFrame(columns=['date','D','R','LV','r','e','l','h','t'])
    eventtype_df = pd.DataFrame(columns=cols)
    eventtype_df.set_index('date',inplace=True)

    # create a new dataframe to do daily counts by analyst
    cols2 = ['date']
    for ak in analyst_hash.keys():
        cols2.append(ak)
    analyst_df = pd.DataFrame(columns=cols2)
    analyst_df.set_index('date',inplace=True)

    for csvfile in list_of_csv_files:

        # read CSV into dataframe
        df = pd.read_csv(csvfile)
        # simplify column names to lowercase and no space
        df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')

        if len(df.mainclass.unique()): # number of mainclasses must be >0

            yyyy = int(csvfile[11:15])
            mm = int(csvfile[15:17])
            start_date = dt.date(yyyy, mm, 1)
            end_date = dt.date(yyyy, mm+1, 1)
            delta = dt.timedelta(days=1)
            while start_date < end_date: # loop over all days in this month
                yyyymmdd = start_date.strftime("%Y%m%d")
                print(yyyymmdd)
                sfiledatetime = pd.to_datetime(df.datetime)
                daycat = df.loc[sfiledatetime.dt.date == start_date] # get rows just for this day
                dayclass_hash = daycat.mainclass.value_counts() # mainclasses for this day
                daysubclass_hash = daycat.subclass.value_counts() # subclasses for this day
                dayanalyst_hash = daycat.analyst.value_counts() # analysts for this day
                # now go through the mainclass and subclass hashes from the overall summary, and make into a new dataframe here
                for ck in class_hash.keys():
                    if ck in dayclass_hash:
                        eventtype_df.at[yyyymmdd, ck] = dayclass_hash[ck] 
                    else:
                        eventtype_df.at[yyyymmdd, ck] = 0
                for sk in subclass_hash.keys():
                    if sk in daysubclass_hash:
                        eventtype_df.at[yyyymmdd, sk] = daysubclass_hash[sk] 
                    else:
                        eventtype_df.at[yyyymmdd, sk] = 0
                for ak in analyst_hash.keys():
                    if ak in dayanalyst_hash:
                        analyst_df.at[yyyymmdd, ak] = dayanalyst_hash[ak] 
                    else:
                        analyst_df.at[yyyymmdd, ak] = 0
                
                start_date += delta # add a day
        #nrows, ncolumns = df.shape
    print(eventtype_df)
    print(analyst_df)
    return eventtype_df, analyst_df

def generate_monthly_wav_csv(mm):
# function generate_monthly_wav_csv(monthdirs)
    print("Reading %s" % mm)
    wavfileslist = sorted(glob.glob(os.path.join(mm, "*")))
    outfile = mm[-13:-8] + 'wavfiles' + mm[-7:-3] + mm[-2:] + '.csv'
    print("...Creating %s" % outfile)
    if os.path.exists(outfile):
        print("CSV file %s already exists. Not overwriting. If you want to overwrite then please delete from the command line" % outfile)
        return
    fptr = open(outfile,'w')
    fptr.write('datetime,duration,wavfilepath,sampling_rate,npts,traceNum,traceID,fixedID,quality_factor,snr,highval,lowval\n') 
    for thiswavfile in wavfileslist:
        if thiswavfile.find('.png')>-1:
            continue
        try:
            st = obspy.read(thiswavfile)
        except:
            print("Cannot read wavfile %s" % thiswavfile )
            continue
        tmp = thiswavfile.split('WAV')
        shortwavfile = 'WAV' + tmp[1] 
        basename = os.path.basename(shortwavfile)
        pos = basename.find('S.')
        #if pos==15: # e.g. 9610-25-1005-41
        #    yyyy = int(basename[0:2]) + 1900
        #elif pos==18:
        #    yyyy = int(basename[0:4])
        #mm=int(basename[pos-13:pos-11])
        #dd=int(basename[pos-10:pos-8])
        #hh=int(basename[pos-7:pos-5])
        #mi=int(basename[pos-5:pos-3])
        #ss=int(basename[pos-2:pos])
        #filetime = dt.datetime(yyyy,mm,dd,hh,mi,ss)
        tracenum = 0
        for tr in st:
            tracetime = tr.stats.starttime
            tr2, quality_factor, snr = qc.compute_metrics(tr)
            duration = tr.stats.npts / tr.stats.sampling_rate
            fptr.write("%s,%s,%8.3f,%3d,%6d,%02d,%s,%s,%d,%.2f,%.2f,%.2f\n" \
                %  (tracetime.datetime, \
                    shortwavfile, \
                    duration, \
                    tr.stats.sampling_rate, tr.stats.npts, \
                    tracenum, tr.id, tr2.id, \
                    quality_factor, snr[0], snr[1], snr[2] ))
            tracenum += 1
            
    fptr.close()

def wav_to_picklefile(wavfile):
    pass

def resp2stationxml(SEISAN_DATA, montserrat=False):
    # see email to obspy users on 2021.07.09
    '''
    We still have the problem that once corrected, an inventory could contain definitions for say MV.MBGH..BHZ and MV.MBGH..HHZ.
    So how to split them into separate files? Or merge everything into one file?
    
    I might also want to add station latitude and longitude from STATION0.HYP too
    '''
    caldir = os.path.join(SEISAN_DATA, 'CAL')
    respfiles = glob(os.path.join(caldir, "RESP.*"))
    for respfile in respfiles:
        xmlfile = respfile.replace('RESP','station') + '.xml'
        if montserrat:
            xmlbase = os.path.basename(xmlfile)
            xmlfile = os.path.join(caldir, xmlfile.replace('MN', 'MV') )
        print('RESP file = ',respfile)
        this_inv = _read_resp(respfile)
        print(this_inv)
        this_inv = inventory_fix_id_mvo(this_inv)
        print(this_inv)
        if os.path.exists(xmlfile): 
            this_inv.write('temp.xml' , format="STATIONXML")
            os.system("cat %s temp.xml > temp2.xml" % xmlfile)
            os.system("rm temp.xml")
            os.system("mv temp2.xml %s" % xmlfile)
        else:
            this_inv.write(xmlfile , format="STATIONXML")


            
def seisanfap2stationxml(fapfile, stationXMLfile):
    # see email to obspy users on 2021.07.09
    pass    

            
def create_event_picklefiles(SEISAN_DATA, DB, YYYY, MM, shortperiod, MAX_FILES_TO_PROCESS=999999):
    """
    Process a year/month of WAV Seisan files into a picklefiles and seismogram plots. This creates
    PICKLE and HTML directories parallel to WAV.
    
    This function came from seisanwavdb2web in PROJECTS/mvocatalog.
    
    To do: should create a separate function just to process on Seisan file at a time? Then this could be
    called as part of a workflow that loads an Sfile, locates the corresponding WAV file, and then processes it to products.
    See wav2picklefile stub above.
    
    Might also want resp2stationxml and seisanfap2stationxml functions. And those could include fixing NSLC ids. See stubs above.
    Then here I can shorten code by just loading stationXML files.
    
    """
    
    # prepare directories. mirroring the seisan/WAV directory, we will store png files in seisan/HTML
    # and pickle files in seisan/PICKLE
    eventdir = os.path.join(SEISAN_DATA, 'WAV', DB, YYYY, MM)
    badseisanfiles = []
    pngdir = eventdir.replace('WAV', 'HTML')
    pickledir = eventdir.replace('WAV', 'PICKLE')
    caldir = os.path.join(SEISAN_DATA, 'CAL')
    if not os.path.exists(pngdir):
        os.makedirs(pngdir)
    if not os.path.exists(pickledir):
        os.makedirs(pickledir)  
        
        
    # Loop over all files found in the seisan/WAV/DB/YYYY/MM directory    
    num_files_processed = 0    
    for root, dirs, files in os.walk(eventdir, topdown=False):
        files = sorted(files)
        for seisanfile in files:

            # we put a limit here for testing purposes
            if num_files_processed >= MAX_FILES_TO_PROCESS:
                print('Have reached MAX_FILES_TO_PROCESS limit')
                break
            
            # check for non-seisan files
            if seisanfile[-4:]=='.png':
                #os.remove(thisfile)
                continue
            if not (seisanfile[0:2]=='19' or seisanfile[0:2]=='20'):
                print('Unrecognized file: %s' % seisanfile)
                continue

            # Okay, we think we have a valid Seisan WAV file to process. But don't load it yet.
            print('Processing %s.' % seisanfile, end=' ')
            rawpngfile = os.path.join(pngdir, seisanfile + '_raw.png')
            correctedpngfile = rawpngfile.replace('_raw', '_corrected')
            picklefile = correctedpngfile.replace('_corrected.png', '.pickle').replace('HTML', 'PICKLE')
            if os.path.exists(picklefile) and os.path.exists(rawpngfile) and os.path.exists(correctedpngfile):
                print('Skipping.')
                continue # files already exist, so skip this event.
            
            # Try to read the Seisan file
            st = Stream()            
            seisanfullpath = os.path.join(root, seisanfile) 
            print('Reading.', end = ' ')
            try:               
                st = read(seisanfullpath)
            except:
                badseisanfiles.append(seisanfile)
                print('ERROR. Could not load.')
                continue # failed to read. log the file to list, and go to next.
            
            if len(st)==0:
                badseisanfiles.append(seisanfile)
                print('ERROR. No traces.')
                continue
            else: 
                print('Success.')

            # raw data - plot all channels
            st.plot(equal_scale=False, outfile=rawpngfile, dpi=100); 

            ######################### START OF CLEAN/CORRECT BLOCK ######################
           
            fix_trace_id(st, shortperiod=shortperiod) 
             
            for tr in st:
                                                   
                ## find inventory
                this_inv = None
                #tr.stats.network = 'MV' 
                #xmlfile = os.path.join(caldir, "station.%s.%s.xml" % (tr.stats.network, tr.stats.station) )
                if tr.stats.channel[1] in 'ES':
                    matchcode = '[ES]'
                elif tr.stats.channel[1] in 'BH':
                    matchcode = '[BH]'
                xmlfiles = glob(os.path.join(caldir, "station.MN.%s..%s*%s.xml" % (tr.stats.station, matchcode, tr.stats.channel[2]) ))
                
                N = len(xmlfiles)
                if N==1:
                    xmlfile = xmlfiles[0]
                    this_inv = read_inventory(xmlfile)   
                    #print('Processing %s with %s' % (tr.id, xmlfile) )
                    process_trace(tr, inv=this_inv)
                else:
                    #print('%d xmlfiles found' % N)
                    process_trace(tr, inv=None)

            
            # if we fixed the NSLCs in the RESP/StationXML file, we would want to fix_trace_id
            # at start of this block, before loop tr in st
            ######################### END OF CLEAN/CORRECT BLOCK ######################
          
            logfile = picklefile.replace('.pickle', '.log')
            with open(logfile,'w') as fout:
                for tr in st:

                    # print trace history
                    pprint(tr.stats, stream=fout)
                    fout.write('\n')
                    """
                    for key in tr.stats:
                        pprint(tr.stats[key], stream=fout)
                        print('\n')
                    """

                    if tr.stats.quality_factor <= 0.0:
                        st.remove(tr)
                    
            if len(st) == 0:
                continue
           
            
            # Write pickle file. This is now the best place to load stream data from in future, so replaces
            # seisan/WAV directory with seisan/PICKLE
            st.write(picklefile, format='PICKLE')
            num_files_processed += 1
            
            # Plot best channels
            # To do: need separate plots for seismic and infrasound channels, as they scale differently
            chosen = choose_best_traces(st, MAX_TRACES=20, include_seismic=True, 
                                        include_infrasound=False, include_uncorrected=False)
            if len(chosen)>0:
                st2 = select_by_index_list(st, chosen)
                st2.plot(equal_scale=True, outfile=correctedpngfile)
                   
    return badseisanfiles


def create_event_spectrograms(SEISAN_DATA, DB, YYYY, MM):                       

    eventdir = os.path.join(SEISAN_DATA, 'PICKLE', DB, YYYY, MM)
    
    for root, dirs, files in os.walk(eventdir, topdown=False):
        files = sorted(files)
        for file in files:
            if not '.pickle' in file:
                continue
                
            print('Processing %s.' % file, end = ' ')
            picklefullpath = os.path.join(root, file) 
            sgramfile = picklefullpath.replace('.pickle', '_sgram.png').replace('PICKLE', 'HTML') 
            sgramfixed = picklefullpath.replace('.pickle', '_sgram_fixed.png').replace('PICKLE', 'HTML')
            if os.path.exists(sgramfile) and os.path.exists(sgramfixed):
                print('Skipping.')
                continue
                     
            st = read(picklefullpath)
      
            iwsobj = IceWeb.icewebSpectrogram(stream=st)
            iwsobj = iwsobj.precompute()
     
            # free scale  
            print('Creating %s.' % sgramfile, end = ' ')            
            titlestr = os.path.basename(sgramfile) 
            chosen = choose_best_traces(st, MAX_TRACES=10)            
            iwsobj.plot(outfile=sgramfile, log=False, equal_scale=False, add_colorbar=True, dbscale=True, title=titlestr, trace_indexes=chosen);

            # fixed scale                        
            titlestr = os.path.basename(sgramfixed)
            clim_in_dB = [-160, -100]
            clim_in_units = [ IceWeb.dB2amp(clim_in_dB[0]),  IceWeb.dB2amp(clim_in_dB[1]) ]
            iwsobj.plot(outfile=sgramfixed, log=False, clim=clim_in_units, add_colorbar=True, dbscale=True, title=titlestr, trace_indexes=chosen);

            # need separate plots for any infrasound channels
            
            # Flatten the spectrogramdata to spectrum, then remove. we use the spectum later.
            print('Computing spectrum.', end = ' ')  
            iwsobj.compute_amplitude_spectrum(compute_bandwidth=True) # adds tr.stats.spectrum
            for tr in iwsobj.stream:
                ampengfft(tr, eventdir)
                tr.stats.pop('spectrogramdata', None) # Remove spectrogramdata as it is too large for picklefile
                #tr.stats.spectrogramdata = {}
            print('Writing enhanced pickle file.')     
            iwsobj.stream.write(picklefullpath, 'PICKLE') # Rewrite pickle file with extra attributes

            
def pickle2csv(SEISAN_DATA, DB, YYYY, MM, MAX_FILES_TO_PROCESS=999999):
    # what metrics do we put into dataframe here. spectraldata. 
    # we want to add any metadata to print on webpages
    # write csvfile replace('.pickle','.csv')                      

    eventdir = os.path.join(SEISAN_DATA, 'PICKLE', DB, YYYY, MM)
    num_files_processed = 0
       
    for root, dirs, files in os.walk(eventdir, topdown=False):
        files = sorted(files)
        for file in files:        
            if num_files_processed >= MAX_FILES_TO_PROCESS:
                print('Have reached MAX_FILES_TO_PROCESS limit.')
                break

            if not '.pickle' in file:
                continue
                
            print('Processing %s.\n- Building metrics dataframe.' % file, end = ' ')
            picklefullpath = os.path.join(root, file) 
            csvfile = picklefullpath.replace('.pickle', '.csv')
            if os.path.exists(csvfile):
                print('Skipping.')
                continue
            
            print('- Building metrics dataframe.', end = ' ')
            st = read(picklefullpath)
            df = pd.DataFrame()
            list_of_rows = []
            for tr in st:
                s = tr.stats
                row = {'id':tr.id, 'starttime':s['starttime'], 
                       'Fs':s.sampling_rate, 'twin':s.duration,
                       'calib':s.calib, 'units':s.units, 
                       'quality':s.quality_factor,
                       'snr':s.snr[0], 'signal':s.snr[1], 'noise':s.snr[2]}
                if 'spectrum' in s: 
                    for item in ['medianF', 'peakF', 'peakA', 'bw_min', 'bw_max']:
                        try:
                            row[item] = s.spectrum[item]
                        except:
                            pass

                if 'metrics' in s:
                    m = s.metrics
                    for item in ['peakamp', 'peaktime', 'energy', 'RSAM_high', 'RSAM_low', 'band_ratio',
                                 'sample_min', 'sample_max', 'sample_mean', 'sample_median', 
                                 'sample_lower_quartile', 'sample_upper_quartile', 'sample_rms', 
                                 'sample_stdev', 'percent_availability', 'num_gaps']:
                                 #'start_gap', 'num_gaps', 'end_gap', 'sum_gaps', 'max_gap', 
                                 #'num_overlaps', 'sum_overlaps', 'num_records', 'record_length', 
                        try:
                            row[item] = m[item]
                        except:
                            pass  
                        
                if 'scipy' in s:
                    row['skewness'] = s.scipy.skewness
                    row['kurtosis'] = s.scipy.kurtosis                        
                        
                list_of_rows.append(row)
            df = pd.DataFrame(list_of_rows)
            df = df.round({'Fs': 2, 'secs': 2, 'quality':2, 'snr':2, 'min':2, 'max':2, 'rms':2})
            
            print('Saving to CSV.')
            
            df.set_index('id')
            df.to_csv(csvfile)
            
            num_files_processed += 1
            
def create_catalog_website(SEISAN_DATA, DB, YYYY, MM):                       

    eventdir = os.path.join(SEISAN_DATA, 'HTML', DB, YYYY, MM)
    allsgramfiles = []
    allseismograms = []
    allwebpages = []
    
    for root, dirs, files in os.walk(eventdir, topdown=False):
        files = sorted(files)
        for file in files:
            if not '_sgram.png' in file:
                continue
                
            print('Processing %s.' % file, end = ' ')
            sgramfullpath = os.path.join(root, file) 
            htmlfile = sgramfullpath.replace('_sgram.png','.html')
            if os.path.exists(htmlfile):
                print('Skipping.')
                continue
            print('Adding to list.')
            allsgramfiles.append(sgramfullpath)                
            allwebpages.append(htmlfile) # need this for previous and next links
            
    for c, sgramfile in enumerate(allsgramfiles):
        seismogramfile = sgramfile.replace('_sgram','_corrected')
        csvfile = sgramfile.replace('_sgram.png','.csv').replace('HTML', 'PICKLE')
        build_event_webpage(allwebpages, c, streamplot=seismogramfile, sgramplot=sgramfile, csv=csvfile)
    

def build_event_webpage(allwebpages, c, streamdf=None, streamplot=None, sgramplot=None, csv=None):
    # put Python code inside % tags, e.g. <% for i in range(10): %>
    
    htmlfile = allwebpages[c]
    print('Creating %s.' % htmlfile)
    lastwebpage = None
    nextwebpage = None
    ind = allwebpages.index(htmlfile)
    if ind>0:
        lastwebpage = allwebpages[ind-1]
    if ind<len(allwebpages)-2:
        nextwebpage = allwebpages[ind+1]
    
    html = """
    <html>
    <head>
    <title>Event webpage</title>
    </head>
    <body>
    <h1>
    """
    html += "%s" % os.path.basename(htmlfile)[0:18]
    
    html += """
    </h1>
    <table border=1>
    <tr>
    """ 
    
    if lastwebpage:
        html += "<td><a href=\"%s\">Previous</a></td>" % os.path.basename(lastwebpage)  
        
        
    html += "<td>"   
    
        
    if sgramplot:
        sgrambase = os.path.basename(sgramplot)
        html += "<h2>Individually scaled</h2><br/><img src=\"%s\"><br/>"% os.path.basename(sgramplot)
        
        sgramplot_fixed = sgramplot.replace('.png','_fixed.png')
        if os.path.exists(sgramplot_fixed):
            html += "<h2>Fixed scale</h2><br/><img src=\"%s\"><br/>"% os.path.basename(sgramplot_fixed)

        
    html += "</td>"         
    
    if nextwebpage:
        html += "<td><a href=\"%s\">Next</a></td>" % os.path.basename(nextwebpage)
        
    html += "</tr>"

     
    if streamplot:
        html += "<tr><td colspan=\"3\"><img src=\"%s\"></td></tr>" % os.path.basename(streamplot)    
     
                
    html += """
    </table>
    """
    
    if csv:
        if os.path.exists(csv):
            df = pd.read_csv(csv)
            df.set_index('id')
            result = df.to_html()
            html += result 
            
    html += """
    </body>
    </html>
    """    

    #print(html)
    
    # Write HTML String to file.html
    fptr = open(htmlfile, "w")
    fptr.write(html)
    fptr.close() 
    
    
def summarize_each_event(SEISAN_DATA, DB, YYYY, MM):
    print("Summarizing WAV files at %s/%s/%s" % (DB, YYYY, MM))
    eventdir = os.path.join(SEISAN_DATA, 'PICKLE', DB, YYYY, MM)
    traceCSVfiles = sorted(glob(os.path.join(eventdir, '[21]*.csv')))
    if not traceCSVfiles:
        return
    wavCSVfile = os.path.join(SEISAN_DATA, 'wav_%s%s%s.csv' % (DB,YYYY,MM))
    summaryLoD = []
    
    # we run this function everytime because we might have added new event files. And it is fast.
    for traceCSVfile in traceCSVfiles: 

        # read pickle file
        picklefile = traceCSVfile.replace('.csv','.pickle')
        print("Processing: %s" % picklefile)
        st = read(picklefile)
        st.filter("highpass", freq=0.5) 
        chosen = choose_best_traces(st, MAX_TRACES=1, include_seismic=True, 
                                    include_infrasound=False, include_uncorrected=False)
        tr = st[chosen[0]]

        # plot best trace
        plt.figure()
        plt.plot(tr.times(), tr.data)
        plt.ylabel(tr.id)   

        # read CSV file
        df = pd.read_csv(traceCSVfile)
        numOfRows = df.shape[0]
        df = df[df["units"] == 'm/s']
        df.sort_values(by=['quality'], inplace=True)
        
        # get median of 10 best rows
        df = df.head(10)
        row = df.median(axis = 0, skipna = True).to_dict()        

        # add columns
        row['WAV']=os.path.basename(traceCSVfile).replace('.csv','')
        row['traces']=numOfRows

        # detect network events
        trig, ontimes, offtimes = detect_network_event(st, sta=0.4, lta=5.0, threshon=6.0, threshoff=0.24)
        print('%s: %d events detected' % (picklefile, len(ontimes)))
        durations = [t['duration'] for t in trig]
        if len(durations)>0:
            bestevent = np.argmax(durations)
            thistrig=trig[int(np.argmax(durations))]
            row['ontime'] = thistrig['time']
            row['offtime']=thistrig['time']+thistrig['duration']            
            for item in ['duration', 'coincidence_sum', 'cft_peak_wmean', 'cft_std_wmean']:
                row[item]=thistrig[item]
            row['detection_quality']=thistrig['coincidence_sum']*thistrig['cft_peak_wmean']*thistrig['cft_std_wmean']

            # add detection lines on plot
            t0 = tr.stats.starttime
            bottom, top = plt.ylim()
            plt.vlines([row['ontime']-t0, row['offtime']-t0], bottom, top )

        plt.show()

        # append row to list of dicts
        #pprint(row)
        summaryLoD.append(row)

        
    summaryDF = pd.DataFrame(summaryLoD)
    summaryDF = summaryDF.rename(columns={'time': 'ontime'})
    summaryDF.drop(summaryDF.filter(regex="Unname"),axis=1, inplace=True)
    print(summaryDF)
    summaryDF = summaryDF.set_index('WAV')
    """
    print(summaryDF)
    for col in summaryDF.columns:
        print('*%s*' % col"""
    summaryDF.to_csv(wavCSVfile, index=True)

"""    
def events_from_channel_detections(st): 
    # event window based on channel detections
    # testeed but never used as is inferior to network event detector
    add_channel_detections(st, lta=5.0, threshon=3.0, threshoff=-0.5)
    for tr in st:
        if 'triggers' in tr.stats:
            for trigpairUTC in tr.stats.triggers:
                duration = trigpairUTC[1]-trigpairUTC[0]
                print(tr.id, trigpairUTC, duration)

    print(' ')            
    [mint, maxt] = get_event_window(st)
    if mint:
        print('Event: ', mint, maxt, maxt-mint)
    else:
        print('No events')    
"""        
def sfiles2csv(SEISAN_DATA, DB, startdate, enddate, csvfile, MAXFILES=99999):
    slist = sorted(get_sfile_list(SEISAN_DATA, DB, startdate, enddate))
    sLoD = []
    for i,sfile in enumerate(slist):
        if i==MAXFILES:
            break
        s = Sfile(sfile, use_mvo_parser=True)
        #s.printEvents()
        d = s.to_dict()
        sLoD.append(d)
    df = pd.DataFrame(sLoD)
    df.to_csv(csvfile)
    
def processREAfiles(SEISAN_DATA, DB, YYYY, MM, MAXFILES=3):
    csvfile=os.path.join(SEISAN_DATA, 'REA', 'rea_%s%s%s.csv' % (DB, YYYY, MM) )
    startdate = dt.datetime(int(YYYY), int(MM), 1)
    if int(MM)<12:
        enddate = dt.datetime(int(YYYY), int(MM)+1, 1)
    else:
        enddate = dt.datetime(int(YYYY+1), 1, 1)
    sfiles2csv(SEISAN_DATA, DB, startdate, enddate, csvfile)
    
def processSeisanYearMonth(SEISAN_DATA, DB, YYYY, MM, MAXFILES=999999):
    
    # Process a month of WAV files. Ultimately this produces a CSV file with 1 line per event.
    processWAVfiles(SEISAN_DATA, DB, YYYY, MM, MAXFILES=MAXFILES) 
    
    # Load the CSV file
    wavCSVfile = os.path.join(SEISAN_DATA, 'wav_%s%s%s.csv' % (DB,YYYY,MM))
    if not os.path.exists(wavCSVfile):
        return
    df = pd.read_csv(wavCSVfile) 
    df['mainclass'] = None
    df['subclass'] = None
    
    # We aim to add a couple of columns from the S-file, and save to this
    reawavCSVfile=os.path.join(SEISAN_DATA, 'reawav_%s%s%s.csv' % (DB, YYYY, MM) )
    
    # Get s-file list
    startdate = dt.datetime(int(YYYY), int(MM), 1)
    if int(MM)<12:
        enddate = dt.datetime(int(YYYY), int(MM)+1, 1)
    else:
        enddate = dt.datetime(int(YYYY+1), 1, 1)
    slist = sorted(get_sfile_list(SEISAN_DATA, DB, startdate, enddate))
    
    for i,sfile in enumerate(slist):
        if i==MAXFILES:
            break
        s = Sfile(sfile, use_mvo_parser=True)
        #s.printEvents()
        d = s.to_dict()
        
        # Are any wavfiles from this S-file in the dataframe?
        for item in ['wavfile1', 'wavfile2']:
            if d[item]:
                wavbase = os.path.basename(d[item])
                ind = df.index[df[item]==wavbase].tolist()
                if ind:
                    df.iloc[ind]['mainclass']=s.mainclass
                    df.iloc[ind]['subclass']=s.subclass
    df.to_csv(csvfile)

def processSeisanYearMonth2(SEISAN_DATA, DB, YYYY, MM, MAXFILES=999999):
    failedWAVfiles=[]
    LoD = []
    
    # We aim to add a couple of columns from the S-file, and save to this
    reawavCSVfile=os.path.join(SEISAN_DATA, 'reawav_%s%s%s.csv' % (DB, YYYY, MM) )
    
    # Get s-file list
    startdate = dt.datetime(int(YYYY), int(MM), 1)
    if int(MM)<12:
        enddate = dt.datetime(int(YYYY), int(MM)+1, 1)
    else:
        enddate = dt.datetime(int(YYYY)+1, 1, 1)
    slist = sorted(get_sfile_list(SEISAN_DATA, DB, startdate, enddate))
    
    for i,sfile in enumerate(slist):
        print('Processing %d of %d: %s' % (i, len(slist), sfile) )
        if i==MAXFILES:
            break
        s = Sfile(sfile, use_mvo_parser=True)
        #s.printEvents()
        d = s.to_dict()
        
        for item in ['wavfile1', 'wavfile2']:
            if d[item]:
                if os.path.exists(d[item]):
                    wavbase = os.path.basename(d[item])
                    if 'MVO' in wavbase:
                        print('Processing ',d[item])
                        eventrow=[]
                        #try:
                        eventrow = processWAV(d[item], shortperiod=False)
                        #except:
                        #    failedWAVfiles.append(wavbase)
                        #else:
                        if eventrow:
                            eventrow['sfile']=os.path.basename(s.path)
                            eventrow['mainclass']=s.mainclass
                            eventrow['subclass']=s.subclass
                            LoD.append(eventrow)
    if LoD:
        df = pd.DataFrame(LoD)
        print('Writing ',reawavCSVfile)
        df.drop(df.filter(regex="Unname"),axis=1, inplace=True)
        df = df.set_index('WAV')       
        df.to_csv(reawavCSVfile, index=True) 
    return failedWAVfiles
    
def processWAVfiles(SEISAN_DATA, DB, YYYY, MM, MAXFILES=3):
    # processWAVfiles('/Users/thompsong/seismo', 'MVOE_', '2005', '05', MAXFILES=10)
    shortperiod = False
    if DB[0:4]=='ASNE':
        shortperiod = True
    #badWAVfiles = create_event_picklefiles(SEISAN_DATA, DB, YYYY, MM, shortperiod, MAX_FILES_TO_PROCESS=MAXFILES)
    create_event_spectrograms(SEISAN_DATA, DB, YYYY, MM)
    pickle2csv(SEISAN_DATA, DB, YYYY, MM, MAX_FILES_TO_PROCESS=MAXFILES)
    create_catalog_website(SEISAN_DATA, DB, YYYY, MM)
    summarize_each_event(SEISAN_DATA, DB, YYYY, MM)
    print(badWAVfiles)
    # could still add ar_pick to this workflow to try to identify P and S waves for regional, VT and hybrid earthquakes
    # this would be at the Trace level
    
def choose_best_events(SEISAN_DATA, NUM_EVENTS, mainclass, subclass):
    # choose_best_events(SEISAN_DATA, 100, 'LV', 't') to choose the best 100 volcano-tectonic events
    pass

def processWAV(wavfile, shortperiod=False, correct_data=False, make_png_files=False):

    WAVDIR = os.path.dirname(wavfile)
    wavbase = os.path.basename(wavfile)
    CALDIR = WAVDIR.replace('WAV', 'CAL')
    
    HTMLDIR = WAVDIR.replace('WAV', 'HTML')
    PICKLEDIR = WAVDIR.replace('WAV', 'PICKLE')     
    if not os.path.exists(HTMLDIR):
        os.makedirs(HTMLDIR)
    if not os.path.exists(PICKLEDIR):
        os.makedirs(PICKLEDIR)  
        
    rawpngfile = os.path.join(HTMLDIR, wavbase + '_raw.png')
    correctedpngfile = rawpngfile.replace('_raw', '_corrected')
    picklefile = os.path.join(PICKLEDIR, wavbase + '.pickle')
    logfile = picklefile.replace('.pickle', '.log')
    sgramfile = os.path.join(HTMLDIR, wavbase + '_sgram.png')
    sgramfixed = sgramfile.replace('_sgram.png', '_sgram_fixed.png')  
    traceCSVfile = picklefile.replace('.pickle', '.csv')
    detectionfile = os.path.join(HTMLDIR, wavbase + '_detection.png')
    wavrow={}
    picklefile_loaded = False
    
    if os.path.exists(picklefile):
        print(picklefile, ' exists')
    else:
        # We need to create the pickle file
        # Try to read the WAV file
        st = Stream()
        print('Processing %s.' % wavbase, end=' ')
        print('Reading.', end = ' ')
        try:               
            st = read(wavfile)
        except:
            print('ERROR. Could not load.')
            return

        if len(st)==0:
            print('ERROR. No traces.')
            return
        else: 
            print('Success.')
        
        # Raw stream plot of all channels
        if not os.path.exists(rawpngfile) and make_png_files:
            st.plot(equal_scale=False, outfile=rawpngfile, dpi=100);

        ######################### START OF CLEAN/CORRECT BLOCK ######################
        print('Cleaning/correcting')
        fix_trace_id(st, shortperiod=shortperiod) 

        for tr in st:

            this_inv = None
            
            if correct_data: # try to find corresponding station XML
                #tr.stats.network = 'MV' 
                #xmlfile = os.path.join(caldir, "station.%s.%s.xml" % (tr.stats.network, tr.stats.station) )
                if tr.stats.channel[1] in 'ES':
                    matchcode = '[ES]'
                elif tr.stats.channel[1] in 'BH':
                    matchcode = '[BH]'
                xmlfiles = glob(os.path.join(CALDIR, "station.MN.%s..%s*%s.xml" % (tr.stats.station, matchcode, tr.stats.channel[2]) ))

                N = len(xmlfiles)
                if N==1:
                    xmlfile = xmlfiles[0]
                    this_inv = read_inventory(xmlfile)   
                    #print('Processing %s with %s' % (tr.id, xmlfile) )
                    process_trace(tr, inv=this_inv)
                else:
                    #print('%d xmlfiles found' % N)
                    process_trace(tr, inv=None)
            else:
                process_trace(tr, inv=None)


        # if we fixed the NSLCs in the RESP/StationXML file, we would want to fix_trace_id
        # at start of this block, before loop tr in st
        ######################### END OF CLEAN/CORRECT BLOCK ######################
            
        # Write pickle file. This is now the best place to load stream data from in future, so replaces
        # seisan/WAV directory with seisan/PICKLE
        print('Writing ',picklefile)
        st.write(picklefile, format='PICKLE')
        picklefile_loaded = True
        
    # save log file
    if not os.path.exists(logfile):
        if not picklefile_loaded:
            st = read(picklefile)
            picklefile_loaded = True
        print('Writing %s' % logfile)
        with open(logfile,'w') as fout:
            for tr in st:
                # print trace history
                pprint(tr.stats, stream=fout)
                fout.write('\n')
                
    # Correct stream plot of best channels
    if not os.path.exists(correctedpngfile) and make_png_files:
        print('Writing ',correctedpngfile)
        if not picklefile_loaded:
            st = read(picklefile)   
        for tr in st:    
            if tr.stats.quality_factor <= 0.0:
                st.remove(tr)               
        chosen = choose_best_traces(st, MAX_TRACES=20, include_seismic=True, 
                                    include_infrasound=False, include_uncorrected=False)
        if len(chosen)>0:
            st2 = select_by_index_list(st, chosen)
            st2.plot(equal_scale=True, outfile=correctedpngfile)                    
        # To do: need separate plots for seismic and infrasound channels, as they scale differently
                
    if not 'energy' in st[0].stats.metrics:
        print('Computing spectrogram data')
        if not picklefile_loaded:
            st = read(picklefile)         
        iwsobj = IceWeb.icewebSpectrogram(stream=st)
        iwsobj = iwsobj.precompute()

        if not os.path.exists(sgramfile) and make_png_files:
            # free scale  
            print('Creating %s.' % sgramfile, end = ' ')            
            titlestr = os.path.basename(sgramfile) 
            chosen = choose_best_traces(st, MAX_TRACES=10)            
            iwsobj.plot(outfile=sgramfile, log=False, equal_scale=False, add_colorbar=True, dbscale=True, title=titlestr, trace_indexes=chosen);

            # fixed scale 
            print('Creating %s.' % sgramfixed, end = ' ')
            titlestr = os.path.basename(sgramfixed)
            clim_in_dB = [-160, -100]
            clim_in_units = [ IceWeb.dB2amp(clim_in_dB[0]),  IceWeb.dB2amp(clim_in_dB[1]) ]
            iwsobj.plot(outfile=sgramfixed, log=False, clim=clim_in_units, add_colorbar=True, dbscale=True, title=titlestr, trace_indexes=chosen);           
            # need separate plots for any infrasound channels

        # Flatten the spectrogramdata to spectrum, then remove. we use the spectrum later.
        print('Computing spectrum.', end = ' ')  
        iwsobj.compute_amplitude_spectrum(compute_bandwidth=True) # adds tr.stats.spectrum
        for tr in iwsobj.stream:
            ampengfft(tr, PICKLEDIR)
            tr.stats.pop('spectrogramdata', None) # Remove spectrogramdata as it is too large for picklefile
            #tr.stats.spectrogramdata = {}
        print('Writing enhanced pickle file.')     
        iwsobj.stream.write(picklefile, 'PICKLE') # Rewrite pickle file with extra attributes

    if os.path.exists(traceCSVfile):
        tracedf = pd.read_csv(traceCSVfile)
    else: 
        print('- Building metrics dataframe.', end = ' ')
        if not picklefile_loaded:
            st = read(picklefile) 
        tracedf = pd.DataFrame()
        list_of_rows = []
        for tr in st:
            s = tr.stats
            row = {'id':tr.id, 'starttime':s['starttime'], 
                   'Fs':s.sampling_rate, 'twin':s.duration,
                   'calib':s.calib, 'units':s.units, 
                   'quality':s.quality_factor,
                   'snr':s.snr[0], 'signal':s.snr[1], 'noise':s.snr[2]}
            if 'spectrum' in s: 
                for item in ['medianF', 'peakF', 'peakA', 'bw_min', 'bw_max']:
                    try:
                        row[item] = s.spectrum[item]
                    except:
                        pass

            if 'metrics' in s:
                m = s.metrics
                for item in ['peakamp', 'peaktime', 'energy', 'RSAM_high', 'RSAM_low', 'band_ratio',
                             'sample_min', 'sample_max', 'sample_mean', 'sample_median', 
                             'sample_lower_quartile', 'sample_upper_quartile', 'sample_rms', 
                             'sample_stdev', 'percent_availability', 'num_gaps']:
                             #'start_gap', 'num_gaps', 'end_gap', 'sum_gaps', 'max_gap', 
                             #'num_overlaps', 'sum_overlaps', 'num_records', 'record_length', 
                    try:
                        row[item] = m[item]
                    except:
                        pass  

            if 'scipy' in s:
                row['skewness'] = s.scipy.skewness
                row['kurtosis'] = s.scipy.kurtosis                        

            list_of_rows.append(row)
        tracedf = pd.DataFrame(list_of_rows)
        tracedf = tracedf.round({'Fs': 2, 'secs': 2, 'quality':2, 'snr':2, 'min':2, 'max':2, 'rms':2})
        print('Saving to CSV.')
        tracedf.set_index('id')
        tracedf.to_csv(traceCSVfile)

        
    # Summarize event
    print('Create a summary row for whole event')
    numOfRows = tracedf.shape[0]
    if correct_data:
        df = tracedf[tracedf["units"] == 'm/s']
    else:
        df = tracedf
    df.sort_values(by=['quality'], inplace=True)
    df = df.head(10) # get median of 10 best rows    
    wavrow = df.median(axis = 0, skipna = True).to_dict()        
    wavrow['WAV']=wavbase
    wavrow['traces']=numOfRows
    
    if detect_event:
        if not picklefile_loaded:
            st = read(picklefile) 
            picklefile_loaded = True
        st.filter("highpass", freq=0.5)    
        trig, ontimes, offtimes = detect_network_event(st, sta=0.4, lta=5.0, threshon=6.0, threshoff=0.24)
        print('%s: %d events detected' % (picklefile, len(ontimes)))
        durations = [t['duration'] for t in trig]
        if len(durations)>0:
            bestevent = np.argmax(durations)
            thistrig=trig[int(np.argmax(durations))]
            wavrow['ontime'] = thistrig['time']
            wavrow['offtime']=thistrig['time']+thistrig['duration']            
            for item in ['duration', 'coincidence_sum', 'cft_peak_wmean', 'cft_std_wmean']:
                wavrow[item]=thistrig[item]
            wavrow['detection_quality']=thistrig['coincidence_sum']*thistrig['cft_peak_wmean']*thistrig['cft_std_wmean']
            
        if make_png_files:
            chosen = choose_best_traces(st, MAX_TRACES=1, include_seismic=True, 
                                include_infrasound=False, include_uncorrected=False)
            tr = st[chosen[0]] 
            plt.figure()
            plt.plot(tr.times(), tr.data)
            plt.ylabel(tr.id)   
            t0 = tr.stats.starttime
            bottom, top = plt.ylim()
            plt.vlines([wavrow['ontime']-t0, wavrow['offtime']-t0], bottom, top )
            plt.savefig(detectionfile)            

    return wavrow

if __name__ == "__main__":
    pass