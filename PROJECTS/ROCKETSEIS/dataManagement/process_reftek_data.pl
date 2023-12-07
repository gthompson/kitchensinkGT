#!/usr/bin/perl
# Workflow for network installed at Beach House on 2018/10/16, for 2018/10/17 launch
# and downloaded 2018/10/17
use File::Basename;

# Step 0: Set some constants
#my $TOPDIR = "/Users/field/Desktop/Glenn/RocketSeis";
my $TOPDIR = "/media/sdd1/ROCKETSEIS";
my $YYYY = "2018";
my $MM = "10";
my $DD = "16";
my $YYYYMMDD = $YYYY.$MM.$DD; # date network was (most recently re-)configured
my $MFILE = "$TOPDIR/CONFIG/network2batch_".$YYYYMMDD.".m";
#my $PFFILE = "$TOPDIR/CONFIG/network".$YYYYMMDD.".pf";
my $PFFILE = "$TOPDIR/CONFIG/network.pf";
my $PARFILE = $PFFILE;
$PARFILE =~ s/pf/par/;
print "topdir = $TOPDIR\nmfile = $MFILE\npf = $PFFILE\npar = $PARFILE\n";
my $TMPDIR = "$TOPDIR/REFTEK_DATA/TMP";

# Step 1: Take each CompactFlash card from RT130 and download with NEO
# This creates a ZIP file of a directory
printf("Have you finished downloading all CF cards with NEO (y/n) ?");
my $answer = <STDIN>;
chomp($answer);
die("Finish downloading CF cards with NEO before proceeding\n") unless ($answer eq 'y' or $answer eq 'Y');

# Step 2: cd to this directory and move the ZIP files
printf("Have you finished moving all of the NEO Zip files to $TOPDIR/REFTEK_DATA/ZIPFILES (y/n) ?");
my $answer = <STDIN>;
chomp($answer);
die("Finish moving all of the NEO Zip files to the ZIPFILES directory before proceeding\n") unless ($answer eq 'y' or $answer eq 'Y');

# Step 3: unzip the ZIP files
chdir("$TOPDIR/REFTEK_DATA/ZIPFILES");
my @allzipfiles = <"*.ZIP">;
chdir($TMPDIR);
if ($#allzipfiles>-1) {
	foreach $myzipfile (@allzipfiles) {
		printf("cp $TOPDIR/REFTEK_DATA/ZIPFILES/$myzipfile $TMPDIR/ \n");
		system("cp $TOPDIR/REFTEK_DATA/ZIPFILES/$myzipfile $TMPDIR/");
		printf("unzip -o $TMPDIR/$myzipfile \n");
		system("unzip -o $TMPDIR/$myzipfile");
		#printf("Hit ENTER to continue\n");
		#<STDIN>;
		system("mv $TOPDIR/REFTEK_DATA/ZIPFILES/$myzipfile $TOPDIR/REFTEK_DATA/ZIPFILES/processed/");
		unlink("$TMPDIR/$myzipfile");
	}
}

# Step 4: Sort raw data - a legacy directory we are no longer adding to
chdir("$TMPDIR/REFTEK_DATA/RAW");
my @DIGITIZER_DIRS = <"$TOPDIR/REFTEK_DATA/RAW/*.*">;
my $this_digitizer_dir;
foreach $this_digitizer_dir (@DIGITIZER_DIRS) {
	my $thisdig = substr($this_digitizer_dir,-4);
	print "$this_digitizer_dir DIGITIZER=$thisdig\n";
	my @daydirs = <"$this_digitizer_dir/*">;
	my $thisdaydir;
	foreach $thisdaydir (@daydirs) {
		my $base = basename $thisdaydir;
		my $targetdir = "$TOPDIR/REFTEK_DATA/SORTED/$base/$thisdig";
		mkdir($targetdir) unless (-e "$TOPDIR/REFTEK_DATA/SORTED/$base");
		#system("mkdir $TOPDIR/REFTEK_DATA/SORTED/$base") unless (-e "$TOPDIR/REFTEK_DATA/SORTED/$base");
		printf("ditto -V $thisdaydir/ $targetdir\n");
		system("ditto -V $thisdaydir/ $targetdir");
		printf("Hit ENTER to continue\n");
		<STDIN>;
	}
}

# Step 5: Sort TMP data - a new directory we are now adding to
chdir($TMPDIR);
my @JDAY_DIRS = <"2??????">;
my $this_jday_dir;
foreach $this_jday_dir (@JDAY_DIRS) {
	print "$this_jday_dir\n";
	my @digitizerdirs = <"$this_jday_dir/????">;
	my $thisdigdir;
	foreach $thisdigdir (@daydirs) {
		#my $base = basename $thisdigdir;
		printf("ditto -V $thisdigdir/ $TOPDIR/REFTEK_DATA/SORTED/$thisdigdir/\n");
		#system("ditto -V $thisdigdir/ $TOPDIR/REFTEK_DATA/SORTED/$thisdigdir/");
	}
}


# Step 6: Create a networkYYYYMMDD.pf file for this experiment/launch
# For example, see $MFILE
# Running this MATLAB M-file creates $PFFILE
# Put this into CONFIG/ directory
#system("matlab -r $MFILE");

# Step 5: Run batch2par on this pf file
system("mv $PARFILE $PARFILE.bak") if (-e $PARFILE);
system("batch2par $PFFILE -m > $PARFILE");
system("sed -i -e 's/rs200spsrs;/1;         /g' $PARFILE");
system("cat $PARFILE");
exit(1);

# Step 7: Process sorted data converting from Reftek format to Miniseed format
# THIS IS THE STEP WHERE WE HAVE TO CHOOSE THE CORRECT PAR FILE
my @daydirs = <"$TOPDIR/REFTEK_DATA/SORTED/20?????">;
my $thisdaydir;
foreach $thisdaydir (@daydirs) {
	print "$thisdaydir\n";

	# Create a list of files to be processed, e.g.
	my $LIST = "tmp_listfiles";
	my $base = basename $thisdaydir;
	system("ls $TOPDIR/REFTEK_DATA/SORTED/$base/????/1/* > $LIST");
	system("ls $TOPDIR/REFTEK_DATA/SORTED/$base/????/9/* >> $LIST");

	# Run rt2ms on these Reftek files #### MODIFY THIS TO CHOOSE CORRECT PARFILE
	system("rt2ms -F $LIST -p $PARFILE -R -L -Y -o $TOPDIR/REFTEK_DATA/MSEED");
	system("rm $LIST");

}



# Step 8: Process Reftek hourly miniseed files into daily Miniseed files
my $mseedyr = "$TOPDIR/REFTEK_DATA/MSEED/Y".$YYYY;
my @daydirs = <$mseedyr."/R???.01">;
my $thisdaydir;
foreach $thisdaydir (@daydirs) {
	print "$thisdaydir\n";

	# Build an Antelope database of these miniseed files
	my $base = basename $thisdaydir;
	my $jday = substr($base,1,3);
	#system("miniseed2db $mseedyr."/$base/*.m $TOPDIR/REFTEK_DATA/DB/db$jday");

	# Make day volumes
	system("miniseed2days -v -U -w $TOPDIR/REFTEK_DATA/DAYFILES/%Y/%j/%{net}.%{sta}.%{loc}.%{chan}.D.%Y.%j $mseedyr/$base/*.m");

}


# Step 9: Build Antelope database out of the day files
my @daydirs = <"$TOPDIR/REFTEK_DATA/DAYFILES/2018/???">;
my $thisdaydir;
foreach $thisdaydir (@daydirs) {
	print "$thisdaydir\n";
	my $base = basename $thisdaydir;
	system("miniseed2db $TOPDIR/REFTEK_DATA/DAYFILES/$YYYY/$base/* REFTEK_DATA/DB/db$base");
}
1;

