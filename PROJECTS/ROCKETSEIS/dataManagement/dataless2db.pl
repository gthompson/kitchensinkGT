#!/usr/bin/perl
foreach $file (glob("response/*.dataless")) {
	$chan = substr($file, -12, -10);
	print "$file = $chan\n";
	if ($chan eq "HH" || $chan eq "BH" || $chan eq "BD" || $chan eq "HD") {
	print "got here\n";
	system("seed2db $file sitedb");
	}
}
