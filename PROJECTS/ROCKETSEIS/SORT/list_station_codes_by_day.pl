#!/usr/bin/perl
use File::Basename qw(basename);
my $daydir = $ARGV[0];
my $basedir = $ARGV[0];
print("SCANNING: $basedir\n");
my @daydirs = glob("$basedir/2??????");
foreach my $daydir (@daydirs) {
	my $day = basename($daydir);
	#print("$day\n");
	my @digitizerdirs = glob("$daydir/????");
	foreach my $digitizerdir (@digitizerdirs) {
		my $digitizer = basename($digitizerdir);
		print("$digitizerdir\n");
		my $pattern = "$daydir/$digitizer/9/*";
		#system("ls $pattern");
		my @filelist = glob("$pattern");
		if ($#filelist == -1) {
			print("-- no 9 files\n");
		}
		foreach my $file (@filelist) {
			$base = basename($file);
			my $result = `strings $file | head -2 | tail -1 | cut -c4-7`;
			chomp($result);
			my $line = 2;
			my $station = $result;
			if (substr($station,0,2) ne "TA" and substr($station,0,2) ne "FI" and  substr($station,0,1) ne "B") {
				$result = `strings $file | head -20 | cut -c4-7`;
				my @parts = split("\n", $result);
				my $found = 0;
				my $line = 0;
				foreach my $part (@parts) {
					$line++;
					if (substr($part,0,2) eq "TA" or substr($part,0,2) eq "FI" or  substr($part,0,1) eq "B") {
						$found = $line;
						$station = $part;
					}
				}
				if (!$found) {	
					print("$file is corrupt?\n");
				}
			} else {
				$station = $result;
			}
			print("$digitizer, $day, $base, $line, $station\n");

		}
	}
}
