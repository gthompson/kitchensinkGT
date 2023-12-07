#!/usr/bin/perl
# Calculate the magnetic anomaly due to a buried magnetized sphere

# USAGE:
$progname = $0;
$progname =~ s/\.pl//g;
die("Usage: $progname DEPTH RADIUS\nOutput file is $progname_DEPTH_RADIUS.csv\n") if ($#ARGV < 1);

$pi = 3.14159;
#sphere is centered at point (0,0,z)
#set depth of center > 0
$z = $ARGV[0]; #meter

#set the radius < $z
$a = $ARGV[1]; #meter

#set magnetization (amp/m)
$minc = 45*$pi/180; #inclination down in rad
$mdec = 0*$pi/180;  #declination east in rad
$mi = 1; #intensity of magnetization (amp/m)
#$mi = 0.25; #intensity of magnetization (amp/m)

#calculate direction cosines for magnetiation
$ml = cos($minc) * cos($mdec);
$mm = cos($minc) * sin($mdec);
$mn = sin($minc);

#calculate the dipole moment assuming a sphere shape
#notice this scales linearly with the
#intensity of magnetization
$dpm = 4/3 * $pi * $a**3 * $mi;

#components of magnetization in x,y,z directions
$mx = $mi * $ml;
$my = $mi * $mm;
$mz = $mi * $mn;

#set earth's field
$einc = 45*$pi/180;  #inclination down in rad
$edec = 0*$pi/180;   #declination east in rad

$el = cos($einc) * cos($edec);
$em = cos($einc) * sin($edec);
$en = sin($einc); 

#proportionality constant is
# magnetic permeability * 1e9 
#to convert of nT
$prop = 400*$pi;
$filename = $progname."_".$z."_".$a.".csv";
print("Saving to $filename\n");
open(FH, '>', $filename) or die $!;
print(FH "x, M\n");
$lineLength = 400; # m
$sampleInterval = 0.1; # m`
$numSamples = $lineLength / $sampleInterval;
for ($j =0; $j<$numSamples; $j++) {

  #px is the northing of the observation point
  #py is the easting of the observation point
  $px = $j*$sampleInterval-$lineLength/2;
  $py = 0;
  $pz = -$z;

  $r2 = $px**2 + $py**2 + $pz**2;
  $r = sqrt($r2);
  $r5 = $r**5;
  $dotprod = $px*$mx + $py*$my + $pz*$mz;

  #components of the anomalous magnetic field
  $bx = $prop*$dpm*(3*$dotprod*$px-$r2*$mx)/$r5;
  $by = $prop*$dpm*(3*$dotprod*$py-$r2*$my)/$r5;
  $bz = $prop*$dpm*(3*$dotprod*$pz-$r2*$mz)/$r5;

  #calculate total anomaly
  #this calculation works for anomaly << total field strength
  $b_total = $el*$bx + $em*$by + $en*$bz;

  #print results 
  print FH "$px, $b_total\n";
}
close(FH);
system("./dipole_plot.py $z $a");
