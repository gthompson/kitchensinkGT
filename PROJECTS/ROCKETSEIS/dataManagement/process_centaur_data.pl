#!/usr/bin/perl
#use Filename::Basename qw(basename);

# Takes data deposited in CENTAUR_DATA/sort and moves it into an SDS directory structure
# DVEL3 data is all in one directory with STAT.NET.LOC.CHAN.YEAR.JDAY.miniseed
# BCHH1 and BCHH2 is already in SDS directory structure and just needs miniseed stripping and then merging with the final directories
# SDS is YEAR/NET/STATION/CHAN.D/NET.STATION.CHAN.D.YEAR.JDAY

# Step 1: Manually Extract data from SD cards into path/to/RocketSeis/CENTAUR_DATA/sort


# Step 2: Process DVEL3 flat structure into SDS
print("STEP 2\n");
#my $CENTAUR_DATA = "/media/sdd1/ROCKETSEIS/CENTAUR_DATA";
my $CENTAUR_DATA = "/raid/data/KennedySpaceCenter/duringPASSCAL/CENTAUR_DATA";
my @years = qw(2018 2019 2020 2021);
my @stations = qw/DVEL3/;
my $net = "FL";
my $thissta;
system("cd $CENTAUR_DATA");
chdir("$CENTAUR_DATA");
foreach $thissta (@stations) {
	my @ddirs = glob("$CENTAUR_DATA/sort/$thissta");
	my $this_ddir;
	foreach $this_ddir (@ddirs) {
		print "$this_ddir\n";
		my $that_dir = $this_ddir;
		$that_dir =~s/sort\///;
		my @miniseeds = glob("$this_ddir/*.miniseed");
		my $thisfile;
		foreach $thisfile (@miniseeds) {
			my @filepathparts = split /\//, $thisfile;
			my ($station,$network,$location,$channel,$year,$jday) = split /\./, @filepathparts[-1];
			if ((substr($channel,0,1) ne "L") and ((substr($channel,1,1) eq "H") or (substr($channel,1,1) eq "D"))) {
				print "@filenameparts\n";
				my $newdir = "$CENTAUR_DATA/sort/$year/$net/$thissta/$channel.D";
				system("mkdir -p $newdir") unless (-e $newdir); 
				my $newfile = join(".", $net, $station, $location, $channel, "D", $year, $jday);
				$newfile =~ s/\.miniseed//;
				$newfile =~ s/sort\///;
				printf("mv $thisfile $newdir/$newfile\n");
				system("mv $thisfile $newdir/$newfile");
				#my $dbname = "$CENTAUR_DATA/db".$year."".$jday;
				#print "miniseed2db $newdir/$newfile $dbname\n";
				#system("miniseed2db $newdir/$newfile $dbname");
			} else {
				printf("unlink($thisfile)\n");
				unlink($thisfile);
			}
		}
	}
}


# Step 3: Move anything from sort directory SDS structures into correct directory & trim miniseed off end
print("STEP 3\n");
my @years = qw(2018 2019 2020 2021);
foreach my $YYYY (@years) {
	my $JDAY_START = "1";
	my $JDAY_END = "366";
	my @stations = qw/BCHH1 BCHH2 DVEL3/;
	my $net = "FL";
	
	my $thissta;
	system("cd $CENTAUR_DATA");
	chdir("$CENTAUR_DATA");
	foreach $thissta (@stations) {
		my @ddirs = glob("$CENTAUR_DATA/sort/$YYYY/$net/$thissta/H??.D");
		my $this_ddir;
		foreach $this_ddir (@ddirs) {
			print "$this_ddir\n";
			my $that_dir = $this_ddir;
			$that_dir =~s/sort\///;
			print "Will move from $this_ddir to $that_dir\n";
			system("mkdir -p $that_dir") unless (-e $that_dir);
			my @miniseeds = glob("$this_ddir/$net*");
			my $thisfile;
			foreach $thisfile (@miniseeds) {
				print "$thisfile\n";
				my $newfile = $thisfile;
				$newfile =~ s/\.miniseed//;
				$newfile =~ s/sort\///;
				$jday = substr($newfile, -3, 3);
				print("mv $thisfile $newfile\n");
				system("mv $thisfile $newfile");
				#print "miniseed2db $newfile db$jday\n";
				#system("miniseed2db $newfile db$jday");
	
			}
		}
	}
}

# Step 4: Shorten the file names for things already in a SDS structure
print("STEP 4\n");
my @years = qw(2018 2019 2020 2021);
foreach my $YYYY (@years) {
	my $JDAY_START = "1";
	my $JDAY_END = "366";
	my @stations = qw/BCHH1 BCHH2 DVEL3/;
	my $net = "FL";
	
	my $thissta;
	system("cd $CENTAUR_DATA");
	chdir("$CENTAUR_DATA");
	foreach $thissta (@stations) {
		my @ddirs = glob("$CENTAUR_DATA/$YYYY/$net/$thissta/H??.D");
		my $this_ddir;
		foreach $this_ddir (@ddirs) {
			print "$this_ddir\n";
			my @miniseeds = glob("$this_ddir/*.miniseed");
			my $thisfile;
			foreach $thisfile (@miniseeds) {
				print "$thisfile\n";
				my $newfile = $thisfile;
				$newfile =~ s/\.miniseed//;
				$jday = substr($newfile, -3, 3);
				print("mv $thisfile $newfile\n");
				system("mv $thisfile $newfile");
				#print "miniseed2db $newfile db$jday\n";
				#system("miniseed2db $newfile db$jday");
	
			}
		}
	}
}

# sds to antelope symlinks
print("STEP 5\n");
my @years = qw(2018 2019 2020 2021);
foreach my $YYYY (@years) {
	my @stations = qw/BCHH1 BCHH2 DVEL3/;
	my $net = "FL";
	
	my $thissta;
	system("cd $CENTAUR_DATA");
	chdir("$CENTAUR_DATA");
	foreach $thissta (@stations) {
		my @ddirs = glob("$CENTAUR_DATA/$YYYY/$net/$thissta/H??.D");
		my $this_ddir;
		foreach $this_ddir (@ddirs) {
			print "$this_ddir\n";
			my @miniseeds = glob("$this_ddir/$net*");
			my $thisfile;
			foreach $thisfile (@miniseeds) {
				my $thisjday = substr($thisfile, -3);
				my $antelopedir = "$CENTAUR_DATA/antelope/$YYYY/$thisjday";
				my @filepathparts = split /\//, $thisfile;
				my $base = $filepathparts[-1];
				my $linkfile = "$antelopedir/$base";
				print("mkdir -p $antelopedir\n");
				system("mkdir -p $antelopedir") unless (-e $antelopedir);
				unless (-e $linkfile) {
					print("ln -s $thisfile $linkfile\n");
					system("ln -s $thisfile $linkfile");
					print "miniseed2db $linkfile $antelopedir/db$jday\n";
					system("miniseed2db $linkfile $antelopedir/db$jday");
				}	
			}
		}
	}
}
