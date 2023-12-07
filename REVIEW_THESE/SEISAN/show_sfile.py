#!/raid/apps/OBSPY/bin/python
# Read S-file (from REA) using Sfile class
# Glenn Thompson 2013/01/06

import os, sys, datetime
sys.path.insert(0, '/home/t/thompsong/Desktop')
import Sfile, mulplt

def main():
    # read in the wav file name from the command line
    sfile = sys.argv[1]
    sf_obj = Sfile.Sfile(sfile)
    dir(sf_obj)
    sf_obj.display()

if __name__ == "__main__":
    main()
