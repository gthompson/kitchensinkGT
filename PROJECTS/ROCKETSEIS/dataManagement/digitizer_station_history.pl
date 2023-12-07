#!/usr/bin/perl
my ($csvfile, $string, $column) = @ARGV;
system("grep $string $csvfile > grepfile.csv");
my $currentstring = "dummy0";
my $newstring = "dummy1";
my $lastgoodline = "";
open(FIN, "<grepfile.csv");
while (my $line=<FIN>) {
	#print($line);
	chomp($line);
	my @parts = split(",",$line);
	if ($#parts==4) {
		$lastgoodline = $line;
		$newstring = $parts[$column-1];
		#print("$newstring ");
		if ($newstring ne $currentstring) {
			print("$line\n");
			$currentstring = $newstring;
		}
	} else {
		#print($#parts." ");
	}
}
close(FIN);
print($lastgoodline."\n");

