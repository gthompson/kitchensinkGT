from libseisGT import add_to_trace_history

def centaur(inputVoltageRange):
    countsPerVolt = 0.4e6 * 40/inputVoltageRange;
    return countsPerVolt

def trillium():
    voltsPerMS = 750; # V / (m/s)
    return voltsPerMS

def infraBSU(HgInThisSensor=0.5):
    # model 0.5" is default
    # 0.1 - 40 Hz flat
    oneInchHg2Pa = 3386.4;
    linearRangeInPa = oneInchHg2Pa * HgInThisSensor;
    selfNoisePa = 5.47e-3;
    voltsPerPa = 46e-6; # from infraBSU quick start guide
    return voltsPerPa

def ChaparralM25():
    # 36 V p2p
    # 0.1 - 200 Hz flat
    selfNoisePa = 3e-3;
    voltsPerPaHighGain = 2.0; # 18 Pa full scale. 
    voltsPerPaLowGain = 0.4; # 90 Pa full scale. recommended for 24-bit digitizers. 
    voltsPerPaVMod = 0.4 * 90/720; # 720 Pa full scale.
    # Volcano mod reduces sensitivity further.
    return voltsPerPaLowGain

countsPerMS = centaur(40.0) * trillium()
countsPerPa40 = centaur(40.0) * infraBSU(0.5)
countsPerPa1 = centaur(1.0) * infraBSU(0.5)
countsPerPaChap = centaur(40.0) * ChaparralM25()

def correctUSFstations(st):
    for tr in st:
        if tr.stats.network=='AM':
            continue
        if tr.stats.sampling_rate>=50:
            #tr.detrend()
            if tr.stats.channel[1]=='H': # a seismic velocity high-gain channel. L for low gain, N for accelerometer
                # trace is in counts. there are 0.3 counts/ (nm/s).
                #tr.data = tr.data / 0.3 * 1e-9 # now in m/s
                calib = countsPerMS
                units = 'm/s'
                                
            if tr.stats.channel[1]=='D': # infraBSU channel?
                #trtr.data = tr.data / countsPerPa40
                # Assumes 0.5" infraBSU sensors at 40V p2p FS
                # But could be 1" or 5", could be 1V p2p FS or could be Chaparral M-25
                units = 'Pa'
                if tr.stats.station=='BCHH1' and tr.stats.channel[2]=='1': # Chaparral M-25
                    calib = countsPerPaChap
                elif tr.stats.station=='BCHH' or tr.stats.station=='BCHH1' or tr.stats.station[0:3]=='SAK' or tr.stats.network=='NU':
                    calib = countsPerPa40
                else:
                    calib = countsPerPa1
            if tr.id=='FL.BCHH4.00.HDF':
                calib=countsPerPaChap
        
            tr.data = tr.data / calib
            tr.stats['calib'] = calib
            tr.stats['units'] = units
            add_to_trace_history(tr, 'calibration_applied')        
if __name__ == "__main__":
    print('calibrations:')
    print('trillium+centaur40 = %f' % countsPerMS)        
    print('infraBSU+centaur40 = %f' % countsPerPa40)        
    print('infraBSU+centaur1 = %f' % countsPerPa1)        
    print('chaparralM25+centaur40 = %f' % countsPerPaChap)        
    print('chaparralM25+centaur1 = %f' % (countsPerPaChap*40))        
