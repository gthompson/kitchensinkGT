#!/usr/bin/perl

# Process the Raspberry Shake data

# Step 1: Get any raspberry shake data
# This assumes that RS3D is plugged into PASSCAL laptop using ethernet cord
# ssh myshake@raspberryshake.local (passwd: shakeme)
# cd /opt/data/archive
# files like 2018/AM/R1E5E/EHZ.D/AM.R1E5E.00.EHZ.D.2018.288
# first find the ip address
if (0) { 
system("ping raspberryshake.local");

## then use scp
system("mv RASPSHAKE_DATA archive");
system("scp -r myshake\@169.254.138.160:/opt/data/archive/ .");
system("mv archive RASPSHAKE_DATA");
}
# this is an SDS directory structure

# Step 2: Build Antelope db
my $YYYY = "2018";
my $JDAY_START = "200";
my $JDAY_END = "290";
my @stations = qw/R1E5E/;
my $net = "AM";

my $thissta;
chdir("RASPSHAKE_DATA");
foreach $thissta (@stations) {
	my @ddirs = <"$YYYY/$net/$thissta/???.D">;
	my $this_ddir;
	foreach $this_ddir (@ddirs) {
		print "$this_ddir\n";
		my @miniseeds = <"$this_ddir/*">;
		my $thisfile;
		foreach $thisfile (@miniseeds) {
			print "$thisfile\n";
			$jday = substr($thisfile, -3, 3);
			#print "miniseed2db $newfile db$jday\n";
			system("miniseed2db $thisfile db$jday");

		}
	}
}

## Step 14: rsync to /Users/gt/shared on my iMac
#rsync -avz . /Volumes/shared/RocketSeis/
