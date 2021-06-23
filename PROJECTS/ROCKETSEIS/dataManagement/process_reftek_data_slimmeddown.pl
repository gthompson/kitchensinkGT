#!/usr/bin/perl
# Workflow for network installed at Beach House on 2018/10/16, for 2018/10/17 launch
# and downloaded 2018/10/17
use File::Basename;

my $TOPDIR = "/raid/data/KennedySpaceCenter/duringPASSCAL";

#print "topdir = $TOPDIR\nmfile = $MFILE\npf = $PFFILE\npar = $PARFILE\n";
my $TMPDIR = "$TOPDIR/REFTEK_DATA/TMP";

chdir($TOPDIR);
# Modify this so it loops over each file in each digitizer directory and checks if 
# an hourly miniseed like MSEED/Y2019/R197.01/2019.197.07.00.00.AB13.1.3.m already 
# exists (for 3rd channel and 07:00 hour of day 2019/197)
# That way, this could convert only new files. Until then, skip this block unless I
# download new data.

# Step 6: Process sorted data converting from Reftek format to Miniseed format
# THIS IS THE STEP WHERE WE HAVE TO CHOOSE THE CORRECT PAR FILE
# We have to create a networkYYYYJJJ.pf file to track each PASSCAL network change
# see ACTIVE_PROJECTS/Rockets/matlab/passcal_create_dbbuild_batchfiles.m
# Running this MATLAB M-file creates a PF FILE for each change date
#system("matlab -r $MFILE");
my @daydirs = sort(<"$TOPDIR/REFTEK_DATA/SORTED/20?????">);
my $thisdaydir;
#my $PFFILE = "$TOPDIR/CONFIG/network.pf";
my $PARFILE;
foreach $thisdaydir (@daydirs) {
	print "$thisdaydir\n";


	# Create a list of files to be processed, e.g.
	my $LIST = "tmp_listfiles";
	my $base = basename $thisdaydir;
	system("ls $TOPDIR/REFTEK_DATA/SORTED/$base/????/1/* > $LIST");
	system("ls $TOPDIR/REFTEK_DATA/SORTED/$base/????/9/* >> $LIST");

	my $PFFILEJDAY = "$TOPDIR/CONFIG/network".$base.".pf";
	if (-e $PFFILEJDAY) {
		my $PARFILEJDAY	= $PFFILEJDAY;
		$PARFILEJDAY =~ s/pf/par/;	
		$PARFILE = $PARFILEJDAY;
		unless (-e $PARFILE) {
			system("mv $PARFILE $PARFILE.bak") if (-e $PARFILE);
			system("batch2par $PFFILEJDAY -m > $PARFILEJDAY");
			system("sed -i -e 's/rs200spsrs;/1;         /g' $PARFILEJDAY");
			system("sed -i -e 's/x1/32/g' $PARFILEJDAY");
		}
		system("cat $PARFILEJDAY");
	}


	# Run rt2ms on these Reftek files #### MODIFY THIS TO CHOOSE CORRECT PARFILE
	print("rt2ms -F $LIST -p $PARFILE -R -L -Y -o $TOPDIR/REFTEK_DATA/MSEED\n");
	system("rt2ms -F $LIST -p $PARFILE -R -L -Y -o $TOPDIR/REFTEK_DATA/MSEED");
	system("rm $LIST");

}


# Step 7: Process Reftek hourly miniseed files into daily Miniseed files
chdir($TOPDIR);
my @yrdirs = sort(<"$TOPDIR/REFTEK_DATA/MSEED/Y20??">);
my $mseedyrdir;
foreach $mseedyrdir (@yrdirs) {
	print("7: $mseedyrdir\n");
	my @daydirs = sort(<"$mseedyrdir/R*01">);
	my $thisdaydir;
	foreach $thisdaydir (@daydirs) {
		print "7: $thisdaydir\n";
		my $base = basename $thisdaydir;

		## Build an Antelope database of these miniseed files

		#my $jday = substr($base,1,3);
		##print("miniseed2db $mseedyr."/$base/*.m $TOPDIR/REFTEK_DATA/DB/db$jday\n");
		###system("miniseed2db $mseedyr."/$base/*.m $TOPDIR/REFTEK_DATA/DB/db$jday");
		my @mseedhourfiles = sort(<"$thisdaydir/*.m">);
		foreach my $mseedhourfile (@mseedhourfiles) {
			print("miniseed2days -v -U -w REFTEK_DATA/DAYFILES/%Y/%j/%{net}.%{sta}.%{loc}.%{chan}.D.%Y.%j $mseedhourfile\n");
			system("miniseed2days -v -U -w REFTEK_DATA/DAYFILES/%Y/%j/%{net}.%{sta}.%{loc}.%{chan}.D.%Y.%j $mseedhourfile");
 		}
	
		# Make day volumes
		#print("miniseed2days -v -U -w REFTEK_DATA/DAYFILES/%Y/%j/%{net}.%{sta}.%{loc}.%{chan}.D.%Y.%j $mseedyrdir/$base/*.m\n");
		#system("miniseed2days -v -U -w REFTEK_DATA/DAYFILES/%Y/%j/%{net}.%{sta}.%{loc}.%{chan}.D.%Y.%j $mseedyrdir/$base/*.m");
		#<STDIN>;	
	}
}

# Step 8: Build Antelope database out of the day files
chdir($TOPDIR);
my @yrdirs = sort(<"$TOPDIR/REFTEK_DATA/DAYFILES/20??">);
my $yrdir;
foreach $yrdir (@yrdirs) {

	# Build an Antelope database of these miniseed files
	my @daydirs = <"$yrdir/???">;
	my $thisdaydir;
	foreach $thisdaydir (@daydirs) {
		print "$thisdaydir\n";
		my $base = basename $thisdaydir;
		my $year = basename $yrdir;
		my $dbname = "db".$year.$base;
		print("miniseed2db $thisdaydir/* $TOPDIR/REFTEK_DATA/DB/$dbname\n");
		system("miniseed2db $thisdaydir/* $TOPDIR/REFTEK_DATA/DB/$dbname");
	}
}
1;

