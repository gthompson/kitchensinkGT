import antelope.datascope as datascope
import matplotlib as mpl
import numpy as np
if 'DISPLAY' in list(os.environ.keys()):
    mpl.use("Agg")
import matplotlib.pyplot as plt
import getopt
import modgiseis

def usage():
    print('Usage: '+sys.argv[0]+' <inputeventdb> <outputimagefile> [<subset_expr> <snum> <enum>]')
    print("""
    where subset_expr is an optional Datascope subset expression 
    Note: your database must have an origin and event table present""")
    print("""\nExample: 
    plot all earthquakes from the AVO catalog with Ml>=1.0 within 20km of Iliamna since 1994/01/01\n
    %s /Seis/Kiska4/picks/Total/Total iliamna_energy.png 'ml >= 1.0 && time > \"1994/01/01\" && deg2km(distance(lat, lon, 60.0319, -153.0918))<20.0'
    """ % (sys.argv[0]))

def main(argv=None):
    try:
            #opts, args = getopt.getopt(argv, 'vh', ['dbpath=', 'outfile=', 'subset_expr=', 'snum=', 'enum='])
            opts, args = getopt.getopt(argv, 'vh')
        if len(args)<2:
            usage()
            sys.exit(2)
        dbpath = args[0]
        outfile = args[1]
        subset_expr = None
        snum = None
        enum = None
        if len(args)>2:
            subset_expr = args[2]
        if len(args)>3:
            snum = args[3]
        if len(args)>4:
            enum = args[4]
    except getopt.GetoptError as e:
            print(e)
            usage()
            sys.exit(2)

    verbose = False

    for o, a in opts:
            if o == "-v":
                    verbose = True
            elif o in ("-h", "--help"):
                    usage()
                    sys.exit()
            else:
                    assert False, "unhandled option"

    if verbose:
        print("dbpath = " + dbpath)
        print("outfile = " + outfile)
        print("subset_expr = " + subset_expr)

    ######## LOAD THE EVENTS 
    # load events from the database dbpath, and apply subset_expr if there is one
    dictorigin = dict();
    dictorigin, numevents = modgiseis.dbgetorigins(dbpath, subset_expr)

    # if we loaded some events, create plots
    if numevents > 0:

        ###### PLOT DATA HERE 

        # Let matplotlib automatically decide where to put date (x-axis) tick marks, and what style of labels to use
        locator = mpl.dates.AutoDateLocator()
        formatter = mpl.dates.AutoDateFormatter(locator)

        # create the figure canvas
        fig1 = plt.figure()

        # add subplot - ml versus time
        ax1 = fig1.add_subplot(311)
        modgiseis.plot_time_ml(ax1, dictorigin, locator, formatter, snum, enum)

        if numevents > 1:

            # Compute bin_edges based on the first and last event times
            bin_edges, snum, enum = modgiseis.compute_bins(dictorigin, snum, enum)
    
            # add subplot - counts versus time
            ax2 = fig1.add_subplot(312)
            modgiseis.plot_counts(ax2, dictorigin, locator, formatter, bin_edges, snum, enum)
        
            # add subplot - energy versus time
            ax3 = fig1.add_subplot(313)
            modgiseis.plot_energy(ax3, dictorigin, locator, formatter, bin_edges, snum, enum)

        ####### SAVE FIGURE

        # save the figure to outfile
        print("- saving to " + outfile)
        fig1.savefig(outfile)

if __name__ == "__main__":
    if len(sys.argv) > 1:
            main(sys.argv[1:])
    else:
        usage()

