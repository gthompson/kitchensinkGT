import os, sys
import numpy as np
from obspy import read, read_inventory
from obspy.io.xseed.core import _read_resp
#from obspy.signal.freqattributes import spectrum
from obspy.imaging.cm import obspy_sequential

sys.path.append('/Users/thompsong/src/volcanoObsPy/LIB')
from libseisGT import eventStatistics
import libMVO

sys.path.append('/Users/thompsong/src/icewebPy')
import IceWeb

            
def seisandb2miniseed(DB, YYYY, MM):
    eventdir = os.path('WAV', DB, YYYY, MM)
    badseisanfiles = []
    pngdir = eventdir.replace('WAV', 'MSEED')
    if not os.path.exists(pngdir):
        os.makedirs(pngdir)
    for root, dirs, files in os.walk(eventdir, topdown=False):
        files = sorted(files)
        for seisanfile in files:
            
            # check for non-seisan files
            if seisanfile[-4:]=='.png':
                #os.remove(thisfile)
                continue
            if not (seisanfile[0:2]=='19' or seisanfile[0:2]=='20'):
                print('Unrecognized file: %s' % seisanfile)
                continue
                
            seisanfullpath = os.path.join(root, seisanfile)    
            rawpngfile = os.path.join(pngdir, seisanfile + '_raw.png')
            correctedpngfile = rawpngfile.replace('_raw', '_corrected')
            rawmseedfile = rawpngfile.replace('.png', '.mseed')
            correctedmseedfile = correctedpngfile.replace('.png', '.mseed')

            
            if not os.path.exists(rawmseedfile):
                try:
                    st = read(seisanfullpath)        
                except:
                    badseisanfiles.append(seisanfile)
            
            if not st:
                print('No traces loaded')
                continue

            # write raw Miniseed file and create plot of Z-components only
            st.write(rawmseedfile)
            stZ = st.select(component='Z')
            stZ.plot(equal_scale=False, outfile=rawpngfile, dpi=100); 

            for tr in st:

                # remove traces with weirdly low sampling rates
                if tr.stats.sampling_rate < 50:
                    st.remove(tr)    
                
                # fix trace ID
                tr.stats.network = 'MV'
                if len(tr.stats.channel==2 and len(tr.stats.location)==1:
                    tr.stats.channel += tr.stats.location
                    tr.stats.location = ''            

                # inventory
                this_inv = None
                respfile = os.path.join('CAL', "RESP.%s" % tr.id)
                xmlfile = os.path.join(pngdir, "inventory.%s.xml" % tr.id)
                #print('StationXML file = ',xmlfile)
                if os.path.exists(xmlfile):
                    this_inv = read_inventory(xmlfile)
                elif os.path.exists(respfile):
                    #print('RESP file = ',respfile)
                    this_inv = _read_resp(respfile)
                    this_inv.write(xmlfile , format="STATIONXML")
                #print(this_inv)
                    
                reconstituteTrace(tr, inv=this_inv)

            # write corrected Miniseed file and create plot of Z-components only
            st.write(correctedmseedfile) 
            stZ = st.select(component='Z')
            stZ.plot(equal_scale=False, outfile=correctedpngfile, dpi=100);
    return badseisanfiles
                       
def process_miniseed_files():                       

    eventdir = os.path('MSEED', DB, YYYY, MM)

    for root, dirs, files in os.walk(eventdir, topdown=False):
        files = sorted(files)
        for file in files:
            if not '_corrected.mseed' in file:
                continue
                
            print('Processing %s' % file)
            mseedfullpath = os.path.join(root, file) 
            st = read(mseedfullpath)
            stZ = st.select(component='Z')    
      
            # spectrogram
            iwsobj = IceWeb.icewebSpectrogram(stream=stZ)
            sgramfile = mseedfullpath.replace('_corrected.mseed', '_sgram.png')
            iwsobj = iwsobj.precompute()
            iwsobj.plot(outfile=sgramfile, log=False, equal_scale=True, add_colorbar=True, dbscale=True)
            # could compute for all channels and then add a component parameter to limit channels plotted
            
            # ampengfft - need to figure best model for this out
            iwsobj.ampengfft() 
            """
            # add an ampengfft method to IceWeb.icewebSpectrogram and modify iwsobj.stream
            # but rather than just replicate old ampengfft which was more like ampengssam,
            # probably better to return full FFT, just as in frequency metrics sampled once-per-minute in GISMO/+iceweb
            # The easy way to do it for events or 10-min or 20-min continuous data chunks is probably:
            1. iwsobj = IceWeb.icewebSpectrogram(st)
            2. iwsobj = iwsobj.precompute()
            3. aef = ampengfft(iwsobj)
            aef here would be a list of dictionaries, one per Trace
            each dict would contain maxamp, energy, ampspec
            do we want to break these down for each minute, for both event and continuous data?
            
            let's also add method for instrument corrected RSAM from my course or the White Island project
            check Miami Lakes too
            Do I also want to archive SSAM?
            """          
            

        if not aef:
            continue
        
            ampengfft(iwsobj)
        
        # write pickle file
        picklefile = mseedfullpath.replace('_corrected.mseed', '.pickle')
        st.write(picklefile, format='PICKLE')
        
        # write an AEF-file
        fptr = open(aeffile, 'w')
        fptr.write('NSLC, amp, eng, f0-1, f1-2, f2-3, f3-4, f4-5, f5-6, f6-7, f7-8, f8-9, f9-10, f10-11, f11-12, f12-13, f13-14, f14-15, f15-16\n')
        for tr in st:
            fptr.write('%s, %4.2e, %4.2e' % (tr.id, tr.stats.peakamp, tr.stats.energy))
            if 'ssam' in tr.stats:
                for x in tr.stats.ssam:
                    fptr.write(', %4.2e' % (x/tr.stats.peakamp))
            fptr.write('\n')
        fptr.close()
        
if __name__ == '__main__':        
    # Montserrat Seisan data processor, one month at a time
    DB = 'MVOE_'
    YYYY = '2005'
    MM = '05'
