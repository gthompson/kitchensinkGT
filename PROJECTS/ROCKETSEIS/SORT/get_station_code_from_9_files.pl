#!/usr/bin/perl
my @filelist = @ARGV;
foreach my $file (@filelist) {
	my $result = `strings $file | head -2 | cut -c4-7`;
	#chomp($result);
	$result =~ s/\R//g;
	$result = substr($result,4,4);
	print("$file: $result\n");
}
