#!/usr/bin/perl

# header
@digitizers = qw/91F8 92B7 9338 9406 98E6 9BC5 9D7C AB13/;
print("yyyyjjj,");
foreach $digitizer (@digitizers) {
	print(" $digitizer,");
}
print("\n");

# rows
my @years = (2018..2021);
for $year (@years) {
	my @days = (1..366); # normal range
	# full limits of PASSCAL RocketSeis experiment from 2018172 to 2021036
	@days = (172..365) if ($year eq 2018);
	@days = (1..36) if ($year eq 2021);

	for $day (@days) {
		next if ($day eq 366 and $year ne 2020); # account for non-leap years
		$day = sprintf("%03d", $day);
		$yd = $year.$day;
		print($yd.",");
		foreach $digitizer (@digitizers) {
			if (-e "$yd/$digitizer") {
				print("    1,");
			} else {
				print("    0,");
			}
		}
		print("\n");
	}
} 
