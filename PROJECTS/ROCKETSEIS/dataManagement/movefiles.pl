#!/usr/bin/perl
@alldirs = glob("???????.????");
foreach $thisdir (@alldirs) {
	print($thisdir."\n");
	chdir($thisdir);
	@secondleveldirs = glob("???????");
	foreach $seconddir (@secondleveldirs) {
		print("\t".$seconddir);
		if (-e "../$seconddir") {
			print(" - Exists\n");
			chdir($seconddir);
			foreach $thirddir (glob("????")) {
				print("\t\t$thirddir");
				if (-e "../../$seconddir/$thirddir") {
					print("- Exists\n");

				} else {
					print(" - Doesnt exist\n");
					print(" - mv $thirddir ../../../organized/$seconddir/$thirddir\n");
					system("mv $thirddir ../../$seconddir/$thirddir");
				}
			}
			chdir("..");
		} else {
			print(" - Doesnt exist\n");
			print(" - mv $seconddir ../../organized/$seconddir\n");
			system("mv $seconddir ../../$seconddir\n");
		}
	}
	chdir("..");
}
