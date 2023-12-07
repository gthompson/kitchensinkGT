#!/usr/bin/perl
die("Usage: $0 YEAR JDAY\n") unless ($#ARGV==1);
my ($year, $jday) = @ARGV;

#$year = 2018;
#$jday = 290;
#$jday = 253;
#$jday = 224;
print("$0 $year $jday\n");
mkdir("EVENTDB");
mkdir("EVENTDB/$year");
my $eventdir = "EVENTDB/$year/$jday";
mkdir($eventdir);
system("cp CENTAUR_DATA/antelope/$year/$jday/*$jday* $eventdir/");
system("cp RASPSHAKE_DATA/$year/AM/*/*.D/*$jday $eventdir/");
system("cp REFTEK_DATA/DAYFILES/$year/$jday/* $eventdir/");
chdir($eventdir);
system("cd $eventdir");
system("miniseed2db * dball".$year.$jday);
1;




