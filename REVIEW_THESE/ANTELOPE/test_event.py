import sys
sys.path.append('F:/GNodePython/bin')
import pylab as plt
import event

n = input('Number of events')
e = event.Event()
#print e.__str__()
e.loadrandom(n)
#print e.__str__()
FontSize = 12
plt.close('all')
plt.figure(1)
plt.axes([0.1, 0.45, 0.5, 0.5])
e.plotlonlat(FontSize)
plt.axes([0.1, 0.25, 0.5, 0.15])
e.plotdepthlon(FontSize)
plt.axes([0.1, 0.05, 0.8, 0.15])    
e.plotdepthtime(FontSize)
plt.axes([0.7, 0.45, 0.2, 0.5])    
e.plotlatdepth(FontSize)
#plt.show()
plt.savefig('F:/deleteme.png')
