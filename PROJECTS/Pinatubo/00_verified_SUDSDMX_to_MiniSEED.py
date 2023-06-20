# %%
import obspy.core as op
import glob
import os
import shutil
import numpy as np
import subprocess
import pandas as pd
import numpy as np
from IPython.display import display, HTML

def read_DMX_file(DMXfile, fix=True, defaultnet=''):
    # DMX read support now (2023) included in ObsPy. Was not available for the Montserrat ASN conversion in 2019.
    # This produces same result as converting DMX to SAC with sud2sac.exe in Win-SUDS, and then reading into ObsPy
    # Has also been tested against sud2gse.exe.
    # sud2msed.exe is messier, because that program completely loses all tr.id info when converting, so all tr.id set to ...
    # Tested on data from Montserrat 1995-6 and Pinatubo 1991
    #
    # ObsPy DMX reader inserts "unk" in place of an unknown network. We do not want this.
    #
    # ObsPy DMX reader reads DMXfile as uint16 and so is all +ve. 
    # sud2sac.exe converts to numbers either side of 0. 
    # Subtracting 2048 from each sample of tr.data corrects data read in using ObsPy DMX reader to match that from SAC
    #
    # Obspy Miniseed writer needs float, not int, so recast as float.
    #
    # Passing fix=False will just run ObsPy DMX reader without applying any corrections.

    print('Reading %s' % DMXfile)
    try:
        st = op.read(DMXfile)
        print('- read okay')
        if fix:
            for tr in st:
                # ObsPy DMX reader sets network to "unk" if blank. We'd rather keep it blank, or 
                # set with explicitly passing defaultnet named argument.
                if tr.stats.network == 'unk':
                    tr.stats.network = defaultnet
                    
                # ObsPy DMX reader falses adds 2048 to each data sample. Remove that here.
                # Also change data type of tr.data from uint to float so we can write trace to MiniSEED later   
                tr.data = tr.data.astype(float) - 2048.0 
    except:
        print('- ObsPy cannot read this demultiplexed SUDS file')        
    return st

def _run_command(cmd, show_terminal_output=True, errorStr=''):
    success = True
    print(cmd)
    result = subprocess.run(cmd, stdout=subprocess.PIPE) #, check=True)  
    #print(result)
    if show_terminal_output:
        if result.stdout:
            #print('\n')
            print("STDOUT:", result.stdout.decode())  # decode the byte-string
            #print('\n')
        if result.stderr:
            print("STDERR:", result.stderr.decode())
            #print('\n')
    if result.stderr:
        return 1, result
    elif errorStr and result.stdout.decode().find(errorStr):
        return 2, result
    else:
        return 0, result
   

def run_IRIG(irigexe, DMXfile, IRIGfile):
    
    if not os.path.exists(DMXfile):
        print("Cannot find %s" % DMXfile)
        return
    
    if IRIGfile==DMXfile:
        print('Refuse to overwrite DMXfile by running irig.exe. IRIGfile must be different to DMXfile.')
        return -1, None
    
    shutil.copyfile(DMXfile, IRIGfile)   
    
    # Use irig.exe to time correct SUDS DMX file
    print('Time correcting ' + IRIGfile)
    cmd = "%s %s" % (irigexe, IRIGfile)

    rtncode, result = _run_command(cmd, errorStr='ERROR')
    
    irig_delta_correction_secs = None
    try:
        index = result.stdout.decode().find('Delta T:')
        irig_delta_correction_secs = float(result.stdout.decode()[index+16:index+27])
    except:
        print('Could not parse IRIG Delta T. Might have a Station not found: IRIG')
    if rtncode:
        print('IRIG failed')
    return rtncode, irig_delta_correction_secs


def SUDSDMX_to_MSED_to_STREAM(sud2msedexe, DMXfile): 
    # SUDS DMX to MSEED using sud2msed. Then ObsPy to read MSEED into Stream
    
    # Create empty Stream object to return
    st = op.Stream() 
    
    if not os.path.exists(DMXfile):
        print("Cannot find %s" % DMXfile)
        return st
    
    # Step 1: Use sud2msed.exe to convert SUDS DMX file to MSEED file
    print('Converting ' + DMXfile + ' to MSEED file')
    MSEEDfile = DMXfile.replace('.DMX','.MSEED')
    cmd = "%s %s %s" % (sud2msedexe, DMXfile, MSEEDfile)
    rtncode, result = _run_command(cmd)
    if rtncode:
        print('Conversion to MSEED failed')
        return st    
    
    # Step 2: Read the MSEED file 
    if not os.path.exists(MSEEDfile):
        print("Cannot find %s" % MSEEDfile)
        return st   

    try:
        print('- Reading ' + MSEEDfile)        
        st = op.read(MSEEDfile)   
    except:
        print('  - FAILED')   
    return st

def SUDSDMX_to_SAC_to_STREAM(sud2sacexe, DMXfile): 
    # SUDS DMX to SAC using sud2sac. Then ObsPy to read SAC into Stream
    # Sacfiles are removed after processing to prevent filling up filesystem.
    
    # Create empty Stream object to return
    st = op.Stream() 
    
    if not os.path.exists(DMXfile):
        print("Cannot find %s" % DMXfile)
        return st
    
    # Step 1: Use sud2sac.exe to convert SUDS DMX file to SAC files
    print('Converting ' + DMXfile + ' to SAC files')
    cmd = "%s %s" % (sud2sacexe, DMXfile)
    rtncode, result = _run_command(cmd)
    #print(rtncode)
    if rtncode:
        print('Conversion to SAC failed')
        return st
    
    # Step 2: Merge the SAC files into a single valid Miniseed file 
    DMXbasename = os.path.basename(DMXfile).replace('.DMX','')
    print('Reading from SAC files')
    
    sacfilelist = glob.glob(DMXbasename + '.sac-???')
    print(sacfilelist)
    if len(sacfilelist) > 0:
        for sacfile in sacfilelist:
            print('- Reading ' + sacfile)
            try:
                sacst = op.read(sacfile);
                #tr.plot();
            except:
                print('  - FAILED')
            else:
                for tr in sacst:
                    tr2 = tr.copy() #.detrend()
                    if not (all(tr2.data==0)): # remove blank channels
                        st = st + tr
    else:
        print('FAILED. No SAC files found')
    return st

def SUDSDMX_to_GSE_to_STREAM(sud2gseexe, DMXfile): 
    # SUDS DMX to GSE using sud2gse. Then ObsPy to read GSE into Stream
    import subprocess
    
    # Create empty Stream object to return
    st = op.Stream() 
    
    if not os.path.exists(DMXfile):
        print("Cannot find %s" % DMXfile)
        return st
    
    # Step 1: Use sud2gse.exe to convert SUDS DMX file to GSEfile
    print('Converting ' + DMXfile + ' to GSE file')
    cmd = "%s %s" % (sud2gseexe, DMXfile)
    rtncode, result = _run_command(cmd)
    if rtncode:
        print('Conversion to GSE failed')
        return st
    
    # Step 2: Read the GSE file 
    GSEfile = DMXfile.replace('.DMX','.GSE')
    if not os.path.exists(GSEfile):
        print("Cannot find %s" % GSEfile)
        return st   

    try:
        print('- Reading ' + GSEfile)        
        st = op.read(GSEfile)   
    except:
        print('  - FAILED')   
    return st
''' 
def checkIfConstant(x):
    bool_equal = np.all(x == x[0])
    value = None
    if bool_equal:
        value=x[0]
    return bool_equal, value
'''

def _compare_Trace_timing(tr0, tr1, errormsglist=[], okay=True):
    if not tr0.stats.npts == tr1.stats.npts:
        errormsglist.append("ERROR: different number of samples: %d %d" % (tr0.stats.npts, tr1.stats.npts) )
        okay = False
    if not tr0.stats.sampling_rate == tr1.stats.sampling_rate:
        errormsglist.append("ERROR: different sampling rates: %.2f Hz %.2f Hz" % 
        (tr0.stats.sampling_rate, tr1.stats.sampling_rate)  ) 
        okay = False
    if not (all(tr0.times()==tr1.times())):
        errormsglist.append("ERROR: different time array")  
        okay = False
        if not (tr0.stats.starttime==tr1.stats.starttime):
            errormsglist.append('ERROR: start times differ by %.3f s' % (tr0.stats.starttime - tr1.stats.starttime) )
        if not (tr0.stats.endtime==tr1.stats.endtime):
            errormsglist.append('ERROR: end times differ by %.3f s' % (tr0.stats.endtime - tr1.stats.endtime) )            
        if not (tr0.times()[-1]==tr1.times()[-1]):
            errormsglist.append('ERROR: different durations %.3f s %.3f s' % (tr0.times()[-1], tr1.times()[-1]) )                    
    #return errormsglist, okay  
    return okay

def _compare_Trace_data(tr0, tr1, errormsglist=[], okay=True):
    if not (all(tr0.data==tr1.data)):
        same=(tr0.data==tr1.data)
        #diff1=tr1.data - tr0.data
        #errormsglist.append('- out of %d samples, %d same' % (tr0.stats.npts, same.sum()) )
        #print('Differences: max: %f, mean: %f' % (np.max(diff1), np.mean(diff1)) )
        if (all(tr0.copy().detrend('linear').data == tr1.copy().detrend('linear').data)):
            errormsglist.append('ERROR: different data array. After detrending, traces are identical')
        else:
            errormsglist.append('ERROR: different data array. After detrending, traces still different')
        okay = False    
    #return errormsglist, okay
    return okay

def _compare_Trace_info(tr0, tr1):
    duration_fraction = np.abs(tr0.times()[-1]-tr1.times()[-1])/tr0.times()[-1]
    thisList = [('id',tr0.id, tr1.id, tr0.id==tr1.id),
              ('npts', tr0.stats.npts, tr1.stats.npts, tr0.stats.npts == tr1.stats.npts),
              ('sampling_rate', tr0.stats.sampling_rate, tr1.stats.sampling_rate, 
               np.abs(tr0.stats.sampling_rate - tr1.stats.sampling_rate)<0.01),
              ('startdate', tr0.stats.starttime.strftime('%Y/%m/%d'), tr1.stats.starttime.strftime('%Y/%m/%d'), 
               tr0.stats.starttime.strftime('%Y/%m/%d') == tr1.stats.starttime.strftime('%Y/%m/%d')),
              ('starttime', tr0.stats.starttime.strftime('%H:%M:%S.%f'), tr1.stats.starttime.strftime('%H:%M:%S.%f'), tr0.stats.starttime == tr1.stats.starttime),
              ('endtime', tr0.stats.endtime.strftime('%H:%M:%S.%f'), tr1.stats.endtime.strftime('%H:%M:%S.%f'), tr0.stats.endtime == tr1.stats.endtime),
              ('duration', tr0.times()[-1], tr1.times()[-1], duration_fraction < 0.001),
              ('data_mean', np.nanmean(tr0.data), np.nanmean(tr1.data), np.nanmean(tr0.data) == np.nanmean(tr1.data) ),
              ('data_std', np.std(tr0.data), np.std(tr1.data), np.std(tr0.data) == np.std(tr1.data) )]
    thisDf = pd.DataFrame(thisList, columns=['metric', 'tr0', 'tr1', 'equal'])
    thisDf.set_index('metric')
    return thisDf

def _summarise_Trace_info(thisDf):
    if thisDf['equal'].all():
        #print('Traces are the same')
        return 0
    else:
        print('Traces are different')
        display(HTML(thisDf.to_html()))
        return 1
    

def _print_errormsglist(errormsglist):
    errormsgset = set(sorted(errormsglist))
    #print(errormsglist)
    for errormsg in list(errormsgset):
        print(errormsg)    

def compare_Streams(st0_in, st1_in):
    # remove any IRIG trace so it does not mess up comparisons
    st0 = st0_in.copy()
    st1 = st1_in.copy()
    remove_IRIG_channel(st0)
    remove_IRIG_channel(st1)
    if not(len(st0)==len(st1)):
        print('ERROR: different number of Trace objects: %d %d' % (len(st0), len(st1)))
    ids0 = [tr.id for tr in st0]
    ids1 = [tr.id for tr in st1]
    
    errormsglist=[]  
    
    ids_match = True
    if not (ids0==ids1):
        print('ERROR: trace ID lists are different')
        ids_match = False
        print(ids0)
        print(ids1) 
        stations0 = sorted([tr.stats.station for tr in st0])
        stations1 = sorted([tr.stats.station for tr in st1])
        if not (stations0==stations1):
            print('ERROR: station lists are different too. Cannot compare trace by trace.\nWill just compare timing info of first trace of each.')
            #errormsglist, okay = _compare_Trace_timing(st0[0], st1[0])  
            #_print_errormsglist(errormsglist)
            thisDf = _compare_Trace_info(st0[0], st1[0])       
            _summarise_Trace_info(thisDf)
            return
        else:
            print('station lists are same')
    else:
        print('trace ID lists are same')
        
    rtncodeSum = 0    
    for i,tr0 in enumerate(st0):
        okay = True
        this_sta = tr0.stats.station
        if ids_match:
            tr1 = st1.select(id=tr0.id)[0]
        else:
            tr1 = st1.select(station=this_sta)[0]
        '''
        okay = _compare_Trace_timing(tr0, tr1, errormsglist=errormsglist, okay=okay)         
        okay = _compare_Trace_timing(tr0, tr1, errormsglist=errormsglist, okay=okay) 
        if okay:
            errormsglist.append('Trace #%d, %s -> %s, is OK' % (i, tr0.id, tr1.id))
        else:
            errormsglist.append('Trace #%d, %s -> %s, FAILED' % (i, tr0.id, tr1.id))
        '''
        thisDf = _compare_Trace_info(tr0, tr1)
        rtncode = _summarise_Trace_info(thisDf)
        rtncodeSum += rtncode
    if rtncodeSum==0:
        print('!!! Stream objects are same !!!')
    #_print_errormsglist(errormsglist)

        
def fix_NSLC_Pinatubo(st, FDSNnet):
    for tr in st:
        sta = tr.stats.station
        tr.stats.network = FDSNnet
        tr.stats.station = sta[:-1]
        tr.stats.channel = 'EH%c' % sta[-1]
        print('Converting %s => %s' % (sta, tr.id))
        
def fix_sampling_rate(st, fs=100.0):
    for tr in st:
        tr.stats.sampling_rate=fs          

def remove_IRIG_channel(st):
    for tr in st:
        if tr.stats.station=='IRIG':
            st.remove(tr) # we do not want to keep the IRIG trace

def remove_blank_traces(st):
    for tr in st:
        if all(tr.data==0):
            st.remove(tr)             

def _get_Stream_stats(st):
    if len(st)>0:
        s = st[0].stats
        return s.starttime, s.sampling_rate
    else:
        return 0, 0
    
def process_one_DMXfile(originalDMXfile, WINSUDSPATH, run_irig=True, run_sud2gse=True, run_sud2sac=True, run_sud2msed=False):
    
    print('********************************************************************')
    print('*** Processing %s *** ' % originalDMXfile)
    print('********************************************************************')
    
    #DMXfile = originalDMXfile
    DMXbasename = os.path.basename(originalDMXfile)
    convertDMXfile = os.path.join(CONVERT_DIR, DMXbasename)
    print('Copying %s to %s' % (originalDMXfile, convertDMXfile))
    shutil.copyfile(originalDMXfile, convertDMXfile) # copy original DMX file to CONVERT dir
    returnDict = {'SUDS_DMX':DMXbasename, 'DMX_starttime':None, 'DMX_fsamp':None}
    #             'irig_rtnCode':None, 'deltaT':None, 'IRIG_starttime':None, 'IRIG_fsamp':None,
    #             'GSE_starttime':None,'GSE_fsamp':None,
    #             'SAC_starttime':None,'SAC_fsamp':None}
    
    # Check that the DMX file headers are readable. If not, return blank Stream and limited resultDF
    badDMX = False
    try:
        st = op.read(originalDMXfile, headonly=True)
        if len(st)==0:
            badDMX = True
    except:
        badDMX = True
    if badDMX:
        print("FAILED to read headers, or zero Trace in Stream. Skipping.")
        return op.Stream(), returnDict
    
    # Get starttime & sampling_rate of original DMX file
    print('\n********* Reading original DMX file *********')
    st_dmx = read_DMX_file(originalDMXfile, fix=True)
    print(st_dmx)
    starttime_dmx, sampling_rate_dmx = _get_Stream_stats(st_dmx)
    #returnHdr.extend(['file','DMX_starttime','DMX_fsamp'])
    #returnRow.extend([DMXbasename, starttime_dmx, sampling_rate_dmx])
    returnDict['DMX_starttime'] = starttime_dmx
    returnDict['DMX_fsamp'] = sampling_rate_dmx

    st_compare = st_dmx
    
    if run_irig:
        
        IRIGfile = os.path.join(CONVERT_DIR, DMXbasename.replace('.DMX','_IRIG.DMX')) # produced by irig.exe
        
        # Apply IRIG.exe to ensure timestamping is correct
        print('\n********* Running IRIG.exe *********')
        returnCode_irig, irig_delta_correction_secs = run_IRIG(os.path.join(WINSUDSPATH, 'irig.exe'), originalDMXfile, IRIGfile)
        
        # Handle return codes / malfunctions of irig.exe
        if returnCode_irig:
            if irig_delta_correction_secs: # None if bad
                print('irig.exe Failed. Trying to use salvaged IRIG Delta T.')
                if np.abs(irig_delta_correction_secs)<3600.0: 
                    st_irig = st_dmx.copy()
                    for tr in st_irig:
                        tr.stats.starttime += irig_delta_correction_secs
                    print(st_irig)
                else:
                    st_irig = op.Stream()
                    print('IRIG Delta T more than 1 hour. Untrustworthy. Skipping this file.')
            else:
                st_irig = op.Stream()
                print('irig.exe Failed. Possibly a sequence error, or no IRIG station. Skipping this file.')
                
        else:  # this means irig.exe probably ran okay
            # Read IRIG-corrected DMX file direct
            convertDMXfile = IRIGfile # this means sud2gse and sud2sac will read IRIG file, not DMX file
            print('\n\n********* Reading IRIG-corrected DMX file directly *********')
            st_irig = read_DMX_file(IRIGfile, fix=True)
            print(st_irig)
            st_compare = st_irig
        starttime_irig, sampling_rate_irig = _get_Stream_stats(st_irig)
        returnDict['irig_rtnCode']=returnCode_irig
        returnDict['deltaT']=irig_delta_correction_secs
        returnDict['IRIG_starttime']=starttime_irig
        returnDict['IRIG_fsamp']=sampling_rate_irig     
    
    # SUDS DMX to GSE using sud2gse. Then ObsPy to read GSE into Stream.
    if run_sud2gse:
        print('\n\n********* Reading via Conversion to GSE *********')        
        st_sud2gse = SUDSDMX_to_GSE_to_STREAM(os.path.join(WINSUDSPATH, 'sud2gse.exe'), convertDMXfile) 
        print(st_sud2gse)
        if len(st_sud2gse)>0 and st_compare:
            #print('read via conversion to GSE complete\n') 
            print('\n*** Comparing DMX reader vs sud2gse')
            compare_Streams(st_compare, st_sud2gse)
        starttime_gse, sampling_rate_gse = _get_Stream_stats(st_sud2gse)
        returnDict['GSE_starttime']=starttime_gse
        returnDict['GSE_fsamp']=sampling_rate_gse         
        
    # SUDS DMX to SAC using sud2sac. Then ObsPy to read GSE into Stream.
    if run_sud2sac:
        print('\n\n********* Reading via Conversion to SAC *********')         
        st_sud2sac = SUDSDMX_to_SAC_to_STREAM(os.path.join(WINSUDSPATH, 'sud2sac.exe'), convertDMXfile) 
        print(st_sud2sac)
        if len(st_sud2sac)>0 and st_compare:
            #print('read via conversion to SAC complete\n') 
            print('\n*** Comparing DMX reader vs sud2sac')
            compare_Streams(st_compare, st_sud2sac)
        starttime_sac, sampling_rate_sac = _get_Stream_stats(st_sud2sac)
        returnDict['SAC_starttime']=starttime_sac
        returnDict['SAC_fsamp']=sampling_rate_sac
        
    # SUDS DMX to MSEED using sud2msed. Then ObsPy to read MSEED into Stream.     
    if run_sud2msed:
        # sud2msed.exe does not work - it hangs or something - and it also loses all tr.id info
        st_sud2msed = SUDSDMX_to_MSED_to_STREAM(os.path.join(WINSUDSPATH, 'sud2msed.exe'), convertDMXfile) 
        print(st_sud2msed)
        if len(st_sud2msed)>0 and st_compare:
            #print('read via conversion to MSEED complete\n') 
            print('\n*** Comparing DMX reader vs sud2msed')
            compare_Streams(st_compare, st_sud2msed) 
        starttime_msd, sampling_rate_msd = _get_Stream_stats(st_sud2msed)
        returnDict['MSD_starttime']=starttime_msd
        returnDict['MSD_fsamp']=sampling_rate_msd
        
    # returns    
    #print(returnHdr, returnRow)
    return st_compare, returnDict

def _Stream_to_SEISANWAV_filename_translation(st, SEISAN_TOP, seisanDBname):
    # Establishing the Seisan WAV filename corresponding to this Stream
    WAVbasename = "%sM.%s_%03d" % (st[0].stats.starttime.strftime('%Y-%m-%d-%H%M-%S'), seisanDBname, len(st))
    WAVfilename = os.path.join(SEISAN_TOP, 'WAV', seisanDBname, 
        st[0].stats.starttime.strftime('%Y'), 
        st[0].stats.starttime.strftime('%m'),
        WAVbasename)    
    return WAVfilename

# Main paths for SUDS DMX
TOPDIR = 'D:\Dropbox\DATA\Pinatubo'
WAVEFORM_DIR = os.path.join(TOPDIR,'WAVEFORMS')
CONVERT_DIR = 'D:\convert'

# paths for files to register into Seisan
os.environ['SEISAN_TOP']='D:\Dropbox\DATA\SEISAN_DB'
SEISAN_TOP = os.getenv('SEISAN_TOP')
seisanDBname = 'PINAT'
#WORpath = os.path.join(os.getenv('SEISAN_TOP'),'WOR')

# Other constants
FDSNnet = 'XB' # assigned by Gale Cox. 1R is code for KSC.
WINSUDSPATH = os.path.join("C:", "Users", "thompsong", "winsuds290", "bin")

# Loop over all files
list_of_returnDict = []
alldirs = glob.glob(os.path.join(WAVEFORM_DIR, '*'))
for thisdir in alldirs:
    allDMXfiles = glob.glob(os.path.join(thisdir, '*.DMX'))
    #allDMXfiles = os.path.join(WAVEFORM_DIR, '199105', '9105011P.DMX')
    for originalDMXfile in allDMXfiles:
        print(originalDMXfile)
        
        # DMX processing
        #st, returnDict = process_one_DMXfile(originalDMXfile, WINSUDSPATH, run_irig=True, run_sud2gse=True, run_sud2sac=True)
        st, returnDict = process_one_DMXfile(originalDMXfile, WINSUDSPATH, run_irig=False, run_sud2gse=False, run_sud2sac=False)
        remove_blank_traces(st)
        if len(st)==0:
            # Bad Stream/event. Nothing to do.
            continue      
        fix_NSLC_Pinatubo(st, FDSNnet)    
        
        # Establishing the Seisan WAV filename corresponding to this SUDS DMX event filename
        WAVfilename = _Stream_to_SEISANWAV_filename_translation(st, SEISAN_TOP, seisanDBname)
        st.write(os.path.join(CONVERT_DIR,os.path.basename(WAVfilename)),'MSEED')
        
        # Registering into Seisan
        #rtncode = register_Stream_into_SeisanDB(st, SEISAN_TOP, seisanDBname) 
        
        # Save WAV filename & return code to dataframe, concatenate, save
        returnDict['SEISAN_WAV_filename']=os.path.basename(WAVfilename)
        #returnDict['registered']=rtncode
        list_of_returnDict.append(returnDict)
        returnDf = pd.DataFrame(list_of_returnDict)
        display(HTML(returnDf.to_html()))
returnDf.to_csv(os.path.join(TOPDIR,'pinatubo_processing.csv')) # Save each time in case of crash before processing all files 

