#!/usr/bin/env python
# coding: utf-8

# # SUDS to MiniSEED wrapper for all Pinatubo waveform data
# This Python3 notebook demonstrates how to convert a PC-SUDS seismic waveform data file into Miniseed
# It requires Win-SUDS Utilities (for programs demux.exe, irig.exe and sud2sac.exe) and ObsPy.

# Do all library imports first
import obspy.core as op
import glob
import matplotlib.pyplot as plt
import os
import shutil

# Set up all paths that will be used
os.chdir('C:\DATA\Pinatubo')
WINSUDSPATH = os.path.join('C:\\', 'winsuds','bin')
#SUDSbasename = 'waveforms/May_1991/9105010W'
WAVEFORM_DIR = 'WAVEFORMS'
CONVERT_DIR = 'convert'
seisanDBname = 'PNTBO'
WORpath = os.path.join(os.getenv('SEISAN_TOP'),'WOR')
demux = os.path.join(WINSUDSPATH, 'demux.exe')
irig = os.path.join(WINSUDSPATH, 'irig.exe')
sud2sac = os.path.join(WINSUDSPATH, 'sud2sac.exe')
sud2msed = os.path.join(WINSUDSPATH, 'sud2msed.exe')
sudsplot = os.path.join(WINSUDSPATH, 'sudsplot.exe')
print('Paths setup okay')

# Functions
def findLatestFile(dirpath):
    lists = os.listdir(dirpath)                                   
    lists.sort(key=lambda fn:os.path.getmtime(os.path.join(dirpath, fn)))
    return os.path.join(dirpath,lists[-1])

def displayFile(dirpath):
    if os.path.exists(dirpath):
        fptr = open(dirpath, 'r')
        str = fptr.read()
        print(str)
    else:
        print(dirpath + ' not found')

# Loop over all files
alldirs = glob.glob(os.path.join(WAVEFORM_DIR, '*'))
for thisdir in alldirs:
    allDMXfiles = glob.glob(os.path.join(thisdir, '*.DMX'))
    for originalWaveformFile in allDMXfiles:

        DMXbasename = os.path.basename(originalWaveformFile)
        DMXfile = os.path.join(CONVERT_DIR, 'original',DMXbasename) # original

        IRIGfile = os.path.join(CONVERT_DIR, 'irig', DMXbasename) # produced by irig.exe
        IRIGcopy = os.path.join(CONVERT_DIR, 'sac', DMXbasename)
        SACbasename = os.path.join(CONVERT_DIR, 'sac', DMXbasename.replace('.DMX', '.sac')) # produced  by sud2sac.exe
        MSEEDfile1 = os.path.join(CONVERT_DIR, 'irig', DMXbasename.replace('.DMX', '.ms')) # produced by sud2msed.exe
        MSEEDfile2 = os.path.join(CONVERT_DIR, 'sac', DMXbasename.replace('.DMX', '.mseed')) # produced by recombining SAC file
        #PHAfile = SUDSbasename + '.PHA' # this might exist if HYPO71 was run to locate the event
        #PUNfile = SUDSbasename + '.PUN' # this might exist if HYPO71 was run and generated a hypocenter

        # Check that orignal DMX file can be read by ObsPy
        shutil.copyfile(originalWaveformFile, DMXfile)
        if not os.path.exists(DMXfile):
            print('FAILED: ' + DMXfile + ' does not exist')
            continue
        try:
            print(DMXfile)
            st = op.read(DMXfile)
            print('- read okay')
            print(st)
            st.plot()
        except:
            print('- ObsPy cannot read this demultiplexed SUDS format')
            continue

        # Time correct the DMX file and convert to SAC files. 
        # Then read in, plot, and attempt top write to MiniSEED.
        # Note: DMX read support now (2023) included in ObsPy. Was not available for the Montserrat ASN conversion in 2019. 
        if not os.path.exists(IRIGfile):
            shutil.copyfile(DMXfile, IRIGfile)
            print('Time correcting ' + IRIGfile)
            os.system(irig + ' ' + IRIGfile)

        print('Reading ' + IRIGfile)
        try:
            st2 = op.read(IRIGfile)
            print('- Success')
        except:
            print('- FAILED')
        else:
            print(st2)

        # Note that we found that sud2msed does not read trace headers, it just 
        # produces one trace from many and also has an incorrect start time. 
        # This is why we use sud2sac, which we have checked against sudsplot.exe
        # and shows the correct trace headers, waveforms and times.
        # Reading IRIG DMX into ObsPy and saving to Miniseed does not work either,
        # as it seems to read duplicate trace IDs, but some without data, and this
        # produces a masked numpy array in tr.data which crashes the MSEED write.
        # SOLUTION: convert to sac first, then read in, append traces, and write.
        print('Converting ' + IRIGfile + ' to SAC files')
        if not os.path.exists(IRIGcopy):
            shutil.copy(IRIGfile, IRIGcopy)
        print(sud2sac + ' '  + IRIGcopy)    
        os.system(sud2sac + ' ' + IRIGcopy) # Convert IRIG DMX file to SAC files with sud2sac.exe 
        print('Reading from SAC files')
        st3 = op.Stream()
        sacfilelist = glob.glob(SACbasename + '-???')
        print(sacfilelist)
        if len(sacfilelist) > 0:
            for sacfile in sacfilelist:
                print('- Reading ' + sacfile)
                try:
                    tr = op.read(sacfile)
                except:
                    print('  - FAILED')
                else:
                    st3 = st3 + tr
            print(st3)
            #st3.plot(equal_scale=False);

            print('Writing SAC file data to %s with ObsPy' % MSEEDfile2)
            try:
                st3.write(MSEEDfile2)
                print('- Success')
            except:
                print('- FAILED')
                continue
        else:
            print('FAILED. No SAC files found matching ' + SACbasename)
            continue

        # Add this Miniseed event waveform file to the Seisan database.
        # The Miniseed file will be moved to [SEISANTOP]/WAV/[SEISANDB]/YYYY/MM
        # A corresponding Seisan S-file (event metadata file) will be created at [SEISANTOP]/REA/[SEISANDB]/YYYY/MM
        # The Seisan programs MAKEREA and AUTOREG are used here. Since they normally require user input, 
        # we create files to simulate this.  
        yy = DMXbasename[0:2]
        mm = DMXbasename[2:4]
        century = '19'
        if yy < '80':
            century = '20'
        yyyy = century + yy
        fptr = open('makerea_wrapper.txt','w')
        fptr.write(seisanDBname + '\n')
        fptr.write(yyyy + mm + '\n')
        fptr.write('\n')
        fptr.write('BOTH\n')
        fptr.close()
        os.system('makerea < makerea_wrapper.txt')

        # Copy the miniseed file to WOR BEFORE running dirf and autoreg
        WORfile = os.path.join(WORpath, os.path.basename(MSEEDfile2))
        cwd=os.getcwd()
        print('cwd = ' + cwd)
        shutil.copyfile(MSEEDfile2, WORfile)

        os.chdir(WORpath)
        os.system('dirf ' + os.path.basename(WORfile))
        fptr = open('autoreg_wrapper.txt','w')
        fptr.write('L\n')
        fptr.write('m\n')
        fptr.write(seisanDBname + '\n')
        fptr.write('gt\n')
        fptr.close()
        os.system('autoreg < autoreg_wrapper.txt')
        os.chdir(cwd)

        # Can now browse event with:
        #    eev yyyymm PNTBO
        # Can also view the latest S-File with the code below
        SFILE = findLatestFile(os.path.join(os.environ['SEISAN_TOP'],'REA',seisanDBname,yyyy,mm))
        print(SFILE)
        displayFile(SFILE)   


# If BUDSPICK generated a phase picks file (file extension *.PHA), add that into the S-file using the Seisan program HYPNOR

def appendPicks(hypnoroutfile, sfilepath):
    if os.path.exists(hypnoroutfile) and os.path.exists(sfilepath):
        pass
    else:
        print('Cannot combine files, one of the inputs does not exist')
    print("Appending " + hypnoroutfile + " to " + sfilepath)
    newlines = list()
    with open(sfilepath, "r") as fs:
        for line in fs:
            #print(line[-2])
            if line[-2] == '7':
                print("- Inserting " + hypnoroutfile)
                #fh = open(hypnoroutfile, 'r')
                with open(hypnoroutfile, 'r') as fh:
                    for hnline in fh:
                        if hnline[-2] != '1':
                            newlines.append(hnline)
                #str = fh.read()
                #fh.close()
                #newlines.append(str)
            else:
                newlines.append(line)
    fs2 = open(sfilepath, "w")
    for newline in newlines:
        fs2.write(newline)
    fs2.close()
    return
            
if os.path.exists(PHAfile):
    fptr = open('hypnor_wrapper.txt','w')
    fptr.write(PHAfile + '\n')
    fptr.write(century + '\n')
    fptr.close()
    displayFile('hypnor_wrapper.txt')
    os.system('hypnor < hypnor_wrapper.txt')
    if os.path.exists('hypnor.out'):
        print('Running HYPNOR on ' + PHAfile)
        appendPicks('hypnor.out', SFILE)
        print(' - Success')
    else:
        print(' - Failed')
else:
    print('No PHA file for ' + SUDSbasename)
displayFile(SFILE)


def HSUMNOR(PUNfile):
    # Based on hsumnor.for in Seisan (which did not work, but this does!)
    # I should add something for the high accuracy lines too, as this only records origin to nearest tenth of second
    outstr = ' ' * 79 + '1'
    agency = 'MVO'
    if os.path.exists(PUNfile):
        #print('0123456789'*8)
        #displayFile(PUNfile)
        #print(' 1996  625 0337 32.9 L  61.588   3.495 15.1  TES 31 1.0 3.2LTES 3.0CTES 3.2LNAO1') # a line from a TEST file in Seisan
        f1 = open(PUNfile, 'r')
        with open(PUNFile, 'r') as f1:
            for line in f1:
                pass
            
        str = line
        #f1.read()
        #f1.close()
        yy = float(str[0:2].strip())
        mm = float(str[2:4].strip())
        dd = float(str[4:6].strip())
        hr = float(str[7:9].strip())
        mi = float(str[9:11].strip())
        sec = float(str[12:17].strip())
        lat = float(str[17:20].strip())
        mlat = float(str[21:26].strip())
        lon = float(str[27:30].strip())
        mlon = float(str[31:36].strip())
        depth = float(str[37:43].strip())
        magstr = str[45:49].strip()
        try:
            magstr = "%3.1f" % float(magstr)
        except:
            magstr = " " * 3
        nphase = float(str[50:53].strip())
        rms = float(str[62:66].strip())
        dlat = float(lat) + float(mlat)/60.0
        dlon = float(lon) + float(mlon)/60.0
        if yy < 80:
            cen = 20
        else:
            cen = 19
        newstr = ' %2d%02d %02d%02d %02d%02d %4.1f L ' % (cen,yy,mm,dd,hr,mi,sec)    
        newstr = newstr + '%7.3f%8.3f%5.1f  %s%3d%4.1f %s' % (dlat, dlon, depth, agency, nphase, rms, magstr)
        outstr = newstr + outstr[len(newstr):]
    return outstr

def insertHYPO71Summary(PUNfile, sfilepath):
    if os.path.exists(PUNfile) and os.path.exists(sfilepath):
        pass
    else:
        print('Cannot combine files, one of the inputs does not exist')
    print("Converting " + PUNfile + " and inserting into " + sfilepath)
    newlines = list()
    hypo71sumstr = HSUMNOR(PUNfile)
    if hypo71sumstr != "":
        newlines.append(hypo71sumstr + "\n")
    with open(sfilepath, "r") as fs:
        for line in fs:
            newlines.append(line)
    fs2 = open(sfilepath, "w")
    for newline in newlines:
        fs2.write(newline)
    fs2.close()
    return

if os.path.exists(PUNfile):
    
    try:
        insertHYPO71Summary(PUNfile, SFILE)
    except:
        print('Could not translate HYPO71 summary file %s' % PUNfile)
        displayFile(PUNfile)
    displayFile(SFILE)
    


def HSUMNOR(PUNfile):
    # Based on hsumnor.for in Seisan (which did not work, but this does!)
    # I should add something for the high accuracy lines too, as this only records origin to nearest tenth of second
    outstr = ' ' * 79 + '1'
    agency = 'MVO'
    if os.path.exists(PUNfile):
        #print('0123456789'*8)
        #displayFile(PUNfile)
        #print(' 1996  625 0337 32.9 L  61.588   3.495 15.1  TES 31 1.0 3.2LTES 3.0CTES 3.2LNAO1') # a line from a TEST file in Seisan
        with open(PUNfile, 'r') as f1:
            for line in f1:
                pass      
        str = line
        yy = float(str[0:2].strip())
        mm = float(str[2:4].strip())
        dd = float(str[4:6].strip())
        hr = float(str[7:9].strip())
        mi = float(str[9:11].strip())
        sec = float(str[12:17].strip())
        lat = float(str[17:20].strip())
        mlat = float(str[21:26].strip())
        lon = float(str[27:30].strip())
        mlon = float(str[31:36].strip())
        depth = float(str[37:43].strip())
        magstr = str[45:49].strip()
        try:
            magstr = "%3.1f" % float(magstr)
        except:
            magstr = " " * 3
        nphase = float(str[50:53].strip())
        rms = float(str[62:66].strip())
        dlat = float(lat) + float(mlat)/60.0
        dlon = float(lon) + float(mlon)/60.0
        if yy < 80:
            cen = 20
        else:
            cen = 19
        newstr = ' %2d%02d %02d%02d %02d%02d %4.1f L ' % (cen,yy,mm,dd,hr,mi,sec)    
        newstr = newstr + '%7.3f%8.3f%5.1f  %s%3d%4.1f %s' % (dlat, dlon, depth, agency, nphase, rms, magstr)
        outstr = newstr + outstr[len(newstr):]
    return outstr

def insertHYPO71Summary(PUNfile, sfilepath):
    if os.path.exists(PUNfile) and os.path.exists(sfilepath):
        pass
    else:
        print('Cannot combine files, one of the inputs does not exist')
    print("Converting " + PUNfile + " and inserting into " + sfilepath)
    newlines = list()
    hypo71sumstr = HSUMNOR(PUNfile)
    if hypo71sumstr != "":
        newlines.append(hypo71sumstr + "\n")
    with open(sfilepath, "r") as fs:
        for line in fs:
            newlines.append(line)
    fs2 = open(sfilepath, "w")
    for newline in newlines:
        fs2.write(newline)
    fs2.close()
    return

if os.path.exists(PUNfile):
    try:
        insertHYPO71Summary(PUNfile, SFILE)
    except:
        print('Could not translate HYPO71 summary file %s' % PUNfile)
        displayFile(PUNfile)
    displayFile(SFILE)
    


# What remains to be done:
# 3. Some DMX files seem to be broken. Are these the ones copied over? Check 9508264I. Check ASNE_/1995/08/ * 1995/08 since these from two different copying processes.
# 6. make the list of number of files each day from the last event label of each day
# 7. Add DOI's for dataset and software & create github project
# 8. process bad???fileList.txt. For bad DMX files, I can check if WVM file exists and if so, convert. Otherwise, see if a valid miniseed file exists in gse_all on newton - but beware that it has not been shifted by 4 hours. Bad PUN files should now work, but I need to figure out how to add them. Would need to write in a way that works even if WVM and DMX file not found.
# 9. Go through whole summary hypo71 list, and try to associate with any event it overlaps with. Easiest to do this in ObsPy. 
# 10. Go through the sheets I am carrying around.
# 
# Modify script so that we do not run makerea each time
# 
# Things for later:
# * do I need to apply 4 hour time shift (+4)?
# * remove duplicate S-files
# * merge with the student teaching cluster database. S-files should have same name (or be within a second). Grab the classification parts and add them. Or go the other way - adding any location/phase pick lines. Beware of 4 hour time difference.
# 
# 

print(DMXfile)


f = os.path.join('fred','bill','1995','07', DMXfile)
os.path.basename(f)[0:-4] + '.mseed'




