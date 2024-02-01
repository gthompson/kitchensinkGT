#!/usr/bin/env python
# coding: utf-8

import os, sys
import numpy as np
import obspy
sys.path.append( os.environ['ANTELOPE'] + '/data/python' )
import antelope.datascope as datascope

def db2Stream(dbpath, stime, etime):
    st = obspy.Stream()
    if os.path.isfile(dbpath):
        dbptr = datascope.dbopen(dbpath, 'r')
    else:
        print(dbpath, ' not found')
        print('pwd=',os.getcwd())
        return
    sepoch = stime.timestamp
    eepoch = etime.timestamp
    #print(stime, etime, sepoch, eepoch)

    dbtableptr = dbptr.lookup(table = 'wfdisc')
    dbtableptr2 = dbptr.lookup(table = 'snetsta')
    dbjoinptr = dbtableptr.join(dbtableptr2)
    nr = dbjoinptr.query("dbRECORD_COUNT")
    #print('Number of records: ',nr)
    st = obspy.Stream()
    # NOT WORKING. MAYBE TRY trload_cssgrp?
    for dbjoinptr.record in range(1):
        #sta, chan, time, samprate, calib = dbtableptr.getv('sta','chan', 'time', 'samprate', 'calib')
        trptr = dbjoinptr.trload_css(f"{sepoch}", f"{eepoch}")
        nrt = trptr.query("dbRECORD_COUNT")
        #print(trptr)
        #trptr.trapply_calib()
        for trptr.record in range(nrt):
            print(trptr.get())
            net, sta, chan, time, samprate, calib = trptr.getv('snet', 'sta','chan', 'time', 'samprate', 'calib')
            #print(trptr.record, sta, chan, time)
            if len(sta)>5:
                sta = sta[0:5]
            #print(trptr.get())
            tr = obspy.Trace(data=np.array(trptr.trdata()))
            tr.stats.delta = 1/samprate
            tr.stats.network = net
            tr.stats.station = sta
            tr.stats.channel = chan[0:3]
            if len(chan)>3:
                tr.stats.location = chan[4:]
            #tr.stats.calib = calib
            tr.stats.starttime = stime
            st.append(tr)
        trptr.free()
    dbptr.free()
    st.merge(method=0, fill_value=0)
    return st

import glob
def mseed2Stream(mseedpath, stime, etime):
    yyyy = stime.strftime('%Y') 
    jday  = stime.strftime('%j')
    files = glob.glob(os.path.join(mseedpath, yyyy, '*', '*', '[HBESC]*.D', f'*{jday}')) # excluding L and other low rate data
    eyyyy = etime.strftime('%Y') 
    ejday  = etime.strftime('%j')
    if ejday != jday:
        files2 = glob.glob(os.path.join(mseedpath, eyyyy, '*', '*', '[HBC]*.D', f'*{ejday}'))   
        files = files + files2
    st = obspy.Stream()
    for f in files:
        this_st = obspy.read(f, format='MSEED')
        this_st.trim(starttime=stime, endtime=etime)
        this_st.merge(method=0, fill_value=0)
        for tr in this_st:
            st.append(tr)
    return st

def get_trlist(st):
    for tr in st:
        tr.stats.network=''
    trlist = sorted([tr.id for tr in st])
    return trlist

def compare_trace_ids(stA, stB, testA=1, testB=2):
    trlistA = get_trlist(stA)
    trlistB = get_trlist(stB)
    same = True
    for trid in trlistA:
        if not trid in trlistB:
            same = False
            print(trid,' in %d, not in %d' % (testA, testB))
    for trid in trlistB:
        if not trid in trlistA:
            same = False
            print(trid,' in %d, not in %d' % (testB, testA))
    if same:
        print('%d and %d are same' % (testA, testB))

def combine_streams(stB, stA):
    appended = False
    for trA in stA:
        found = False
        for trB in stB:
            if trA.stats.station == trB.stats.station and trA.stats.location == trB.stats.location and trA.stats.channel == trB.stats.channel:
                if trA.stats.starttime >= trB.stats.startttime and trA.stats.endtime <= trB.stats.endtime:
                    found = True
                    break
        if not found:
            stB.append(trA)
            appended = True
    if appended:
        stB.merge(method=0, fill_value=0)
               

if __name__ == "__main__":
    import time

    # Test 1
    launchtime = obspy.UTCDateTime(2022,11,6,13,38,0)
    pretrig = 60
    posttrig = 200
    stime = launchtime - pretrig
    etime = launchtime + posttrig
    time1 = time.time()
    dbpath = f'/home/thompsong/work/PROJECTS/KSC_EROSION/db/db{stime.strftime("%Y%m%d")}'
    st1 = db2Stream(dbpath, stime, etime)
    time2 = time.time()
    print('db2Stream took ',time2-time1)
    print(st1)
    #st1.write('test.mseed','MSEED')
    pngfile = 'test1_' + os.path.basename(sys.argv[0]).replace('.py', '.png')
    print('Writing ',pngfile)
    #st1.plot(equal_scale=False, outfile=pngfile);

    # Test 2
    time1 = time.time()
    mseedpath = '/home/thompsong/work/PROJECTS/KSC_EROSION/SDS'
    st2 = mseed2Stream(mseedpath, stime, etime)
    time2 = time.time()
    print('mseed2Stream took ',time2-time1)
    pngfile = 'test2_' + os.path.basename(sys.argv[0]).replace('.py', '.png')
    print('Writing ',pngfile)
    #st2.plot(equal_scale=False, outfile=pngfile);
    print(st2)

    # Test 3
    from obspy.clients.filesystem.sds import Client
    sdsclient = Client(mseedpath)
    #sdsclient.get_all_nslc(sds_type=None, datetime=stime)
    time1 = time.time()
    st3 = sdsclient.get_waveforms("*", "*", "*", "[HD]*", stime, etime)
    time2 = time.time()
    print('sdsclient.get_waveforms took ',time2-time1)
    pngfile = 'test3_' + os.path.basename(sys.argv[0]).replace('.py', '.png')
    print('Writing ',pngfile)
    #st3.plot(equal_scale=False, outfile=pngfile); 
    print(st3)

    # Test 4
    import header
    paths = header.setup_environment()
    import SDS
    paths['SDS_TOP'] = os.path.join(paths['outdir'], 'SDS')
    thisSDSobj = SDS.SDSobj(paths['SDS_TOP'])
    time1 = time.time()
    thisSDSobj.read(stime, etime, speed=1)
    time2 = time.time()
    print('SDS class took ',time2-time1) 
    pngfile = 'test4_' + os.path.basename(sys.argv[0]).replace('.py', '.png')
    print('Writing ',pngfile) 
    st4 = thisSDSobj.stream
    #st4.plot(equal_scale=False, outfile=pngfile);
    print(st4)
            
    compare_trace_ids(st1, st2, testA=1, testB=2)
    compare_trace_ids(st3, st2, testA=3, testB=2)
    compare_trace_ids(st4, st2, testA=4, testB=2)