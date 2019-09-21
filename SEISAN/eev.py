#!/raid/apps/OBSPY/bin/python
# Make a list of Sfiles or wavfiles
# Glenn Thompson 2013/01/06
# Usage: eev filenr.lis

import sys
sys.path.insert(0, '/home/t/thompsong/Desktop') # path to my Python codes?
import Seisan_Catalog, StreamGT
from os import walk

def main():
    # read in the wav file name from the command line
    filenrlis = sys.argv[1]
    files = []
    file=open(filenrlis,'r') 
        
    lines = file.readlines()

    linenum = 0 
    while True:
        thisfile = line[6:].strip()
        choice = raw_input(line + "? ")
        if choice == 'n':
            continue
    
        elif choice == 's':
            show_sfile(thisfile)

        elif choice == 'p':
            mulplt(thisfile)

        elif choice == 'h':
            print "Menu:"
            print "\tf = first event"
            print "\tn = next event"
            print "\tp = previous event"
            print "\tl = last event"
            print "\tw = plot waveform file(s)"
            print "\ts = show Sfile"
            print "\ta = show AEFfile"
            print "\tq = quit"
	    print " "

        elif choice == 'q':
            break


        
