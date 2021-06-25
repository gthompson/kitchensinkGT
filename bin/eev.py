#!/usr/bin/env python
# Make a list of Sfiles or wavfiles
# Glenn Thompson 2013/01/06
# Usage: eev filenr.lis

import sys
sys.path.insert(0, '/home/t/thompsong/Desktop') # path to my Python codes?
import Seisan_Catalog, StreamGT
from os import walk

def main():
    filenrlis = sys.argv[1]
    file = open(filenrlis,'r') 
        
    lines = file.readlines()
    nlines = len(lines)

    linenum = 0 
    while True:
        thisfile = lines[linenum][6:].strip()
        choice = input(line + "? ")

        if choice == 'n':
            linenum += 1
            if linenum>=nlines:
                linenum=nlines-1
    
        elif choice == 'b':
            linenum -= 1
            if linenum<0:
                linenum=0

        elif choice == 's':
            Seisan_Catalog.show_sfile(thisfile)

        elif choice == 'p':
            StreamGT.mulplt(thisfile)

        elif choice == 'h':
            print("Menu:")
            print("\tf = first event")
            print("\tn = next event")
            print("\tb = previous event")
            print("\tl = last event")
            print("\tw = plot waveform file(s)")
            print("\ts = show Sfile")
            print("\ta = show AEFfile")
            print("\tq = quit")
        print(" ")

        elif choice == 'q':
            break


        
