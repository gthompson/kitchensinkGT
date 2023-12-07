import os, sys
HOME = os.getenv('HOME')
sys.path.append('%s/src/icewebPy' % HOME)
import IceWeb
os.getcwd()
datadir = os.path.join(os.getenv('HOME'), 'DATA', 'MVO', 'PNG')
seisanfile = '%s/MVOE_/2005/05/2005-05-01-0149-37S.MVO___025.mseed' % datadir
st = read(seisanfile)
stZ = st.select(component='Z')
pngfile = os.path.basename(seisanfile) + '.png'
obj = IceWeb.icewebSpectrogram(stream=stZ)
print(obj.t_range())
obj.precompute()
sgramfile = os.path.basename(seisanfile) + '_sgram.png'
obj.plot(outfile=sgramfile, log=False, \
         equal_scale=False, add_colorbar=True, title="Unscaled")        
obj.plot(outfile=sgramfile.replace('sgram','sgram_scaled'), log=False, \
         equal_scale=True, add_colorbar=True, title="Equally Scaled")
obj.plot(outfile=sgramfile.replace('sgram','sgram_fixed'), log=True, \
         clim=[1e-8, 1e-5], add_colorbar=True, title="Fixed", fmin=0.1, fmax=40.0, dbscale=True, cmap=obspy_sequential)
            
