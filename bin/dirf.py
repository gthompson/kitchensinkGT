#!/raid/apps/OBSPY/bin/python
# Make a list of Sfiles or wavfiles
# Glenn Thompson 2013/01/06
# Usage: dirf.py parent-directory filepattern  > filenr.lis
# e.g. dirf.py /raid/data/suds/DMX/1995/07 gse
import sys
from os import walk
#sys.path.insert(0, '/home/t/thompsong/Desktop')

def main():
    # read in the wav file name from the command line
    fileroot = sys.argv[1]
    files = []
    for (dirpath, dirnames, filenames) in walk(fileroot):
            files.extend(filenames)
            break

    i = 1
    for f in sorted(files):
        if sys.argv[2] in f:
            fullpath = fileroot + "/" + f
            print("#%06d\t%s" % (i, fullpath))
            i += 1
    
if __name__ == "__main__":
    main()
