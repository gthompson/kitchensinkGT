#!/usr/bin/env python
import obspy
import os, sys

def usage():
    print("Usage: python %s seismicfile" % sys.argv[0])

def quicklook(seismicfile):
    if os.path.isfile(seismicfile):
        try:
            st = obspy.read(seismicfile)
            print(st)
            print(st[0].stats)
            st.plot(equal_scale=False)
        except:
            print("Cannot read: %s" % seismicfile)
    else:
        print("File not found: %s" % seismicfile)

if len(sys.argv)==2:
    seismicfile = sys.argv[1]
    quicklook(seismicfile)
elif len(sys.argv)==1:
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        seismicfile = filedialog.askopenfilename()
        quicklook(seismicfile)
    except:
        #print('tkinter failed')
        usage()
else:
    usage()

