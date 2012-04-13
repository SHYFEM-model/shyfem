#!/usr/bin/perl
#
# reads custom file and converts points to GRD format
#
# You will have to handle the format of the input file.
# Here its is assumed that each line contains "x y z" values.
# If the units of (x,y) are not in meters you will have to adjust $xyfact.
# If the units of z are not in meters you will have to adjust $zfact.
# For example, if your depth values are negative, use $zfact=-1.
#
#-----------------------------------------------------------------

$xyfact = 1.;		# scale x/y values
$zfact = 1.;		# scale depth values

#-----------------------------------------------------------------

$type = 2;

while(<>) {

  chomp;
  next if /^\s*$/;		# empty line
  s/,/ /g;			# change , to space
  s/;/ /g;			# change ; to space

  ($x,$y,$z) = split;		# handle each line of input file

  $n++;
  $x *= $xyfact;
  $y *= $xyfact;
  $z *= $zfact;

  print "1 $n $type $x $y $z\n";

}

#-----------------------------------------------------------------

