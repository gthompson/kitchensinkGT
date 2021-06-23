#!/usr/bin/perl
#
@zipfiles = glob("*.ZIP");
foreach $zipfile (@zipfiles) {
	system("unzip -B $zipfile");
	system("mv $zipfile processed/");
}

