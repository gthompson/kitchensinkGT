from obspy.clients.nrl import NRL
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.core import UTCDateTime

from libseisGT import add_to_trace_history

def centaur(inputVoltageRange):
    countsPerVolt = 0.4e6 * 40/inputVoltageRange;
    return countsPerVolt

def trillium():
    voltsPerMS = 754; # V / (m/s)
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
            tr.stats['sensitivity'] = calib
            tr.stats['units'] = units
            add_to_trace_history(tr, 'calibration_applied')        


def make_inv(net, sta, loc, chans, datalogger='Centaur', sensor='TCP', Vpp=40, fsamp=100, lat=0.0, lon=0.0, elev=0.0, depth=0.0, sitename='', ondate=UTCDateTime(1970,1,1), offdate=UTCDateTime(2025,12,31)):
    nrl = NRL('http://ds.iris.edu/NRL/')
    if datalogger == 'Centaur':
        if Vpp==40:
            datalogger_keys = ['Nanometrics', 'Centaur', '40 Vpp (1)', 'Off', 'Linear phase', "%d" % fsamp]
        elif Vpp==1:
            datalogger_keys = ['Nanometrics', 'Centaur', '1 Vpp (40)', 'Off', 'Linear phase', "%d" % fsamp]
    elif datalogger == 'RT130':
        datalogger_keys = ['REF TEK', 'RT 130 & 130-SMA', '1', "%d" % fsamp]
    else:
        print(datalogger, ' not recognized')
        print(nrl.dataloggers[datalogger])
    print(datalogger_keys)

    if sensor == 'TCP':
        sensor_keys = ['Nanometrics', 'Trillium Compact 120 (Vault, Posthole, OBS)', '754 V/m/s']
    elif sensor == 'L-22':
        sensor_keys = ['Sercel/Mark Products','L-22D','2200 Ohms','10854 Ohms']
    elif sensor == 'Chap':
        sensor_keys = ['Chaparral Physics', '25', 'Low: 0.4 V/Pa']
    else:
        print(sensor, ' not recognized')
        print(nrl.sensors[sensor])
    print(sensor_keys)



    response = nrl.get_response(sensor_keys=sensor_keys, datalogger_keys=datalogger_keys)
    print(response)
    print(response.instrument_sensitivity)
    channels = []
    for chan in chans:
        channel = Channel(code=chan,
                      location_code=loc,
                      latitude=lat,
                      longitude=lon,
                      elevation=elev,
                      depth=depth,
                      sample_rate=fsamp,
                      start_date=ondate,
                      end_date=offdate,
                      )
        channel.response = response
        channels.append(channel)
    station = Station(code=sta,
                      latitude=lat,
                      longitude=lon,
                      elevation=elev,
                      creation_date=ondate,
                      site=Site(name=sitename),
                      channels=channels,
                      start_date=ondate,
                      end_date=offdate,
                      )

    network = Network(code=net,
                     stations=[station])
    inventory = Inventory(networks=[network], source="demo")
    return inventory




if __name__ == '__main__':
    inv = make_inv('FL', 'BCHH', '', ['HHZ', 'HHN', 'HHE'], datalogger='Centaur', sensor='TCP', Vpp=40, fsamp=100, sitename='Beach House original', ondate=UTCDateTime(2016,2,24) )
    #inv.write("BCHH_original_seismic.sml", format="stationxml", validate=True)
    inv = make_inv('FL', 'BCHH', '', ['HDF'], datalogger='Centaur', sensor='Chap', Vpp=40, fsamp=100, sitename='Beach House Sonic', ondate=UTCDateTime(2017,8,1) )
    #inv.write("BCHH_sonic_Chap.sml", format="stationxml", validate=True) # Replace PA with Pa and 400,000 in InstrumentSensivitity Value with 160,000
   
    print('calibrations:')
    print('trillium+centaur40 = %f' % countsPerMS)        
    print('infraBSU+centaur40 = %f' % countsPerPa40)        
    print('infraBSU+centaur1 = %f' % countsPerPa1)        
    print('chaparralM25+centaur40 = %f' % countsPerPaChap)        
    print('chaparralM25+centaur1 = %f' % (countsPerPaChap*40))        
