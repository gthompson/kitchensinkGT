import obspy as op
import matplotlib.pyplot as plt
st=op.read('/media/data/KennedySpaceCenter/wfdata/rocketevents/20160322_ULA_Atlas_V/BCHH_0383_20160323_0305.seed')


st.plot(outfile='/var/www/usfseismiclab.org/rocket_seismograms.png')
#st.spectrogram()


fp = open('/var/www/usfseismiclab.org/html/cesil2.html','w+')
fp.write(
    '''<html>
    <head>
    <title>Rocket event webpage</title>
    </head>
    <body>
    <h1>rocket event on 2016/03/22</h1>
    <img src="rocket_seismograms.png">
    </body>
    </html>
    ''')
fp.close()
