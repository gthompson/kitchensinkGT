import sys
sys.path.append('F:/GNodePython/bin')
import pylab as plt
import event

n = input('Number of events')
e = event.Event()
#print e.__str__()
e.loadrandom(n)
#print e.__str__()

e.writetoxml()
