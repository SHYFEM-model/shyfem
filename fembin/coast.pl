#!/usr/bin/perl -s
#
# reads custom file and converts points to coast line
#
# You will have to handle the format of the input file.
# Here its is assumed that each line contains "x y" values.
# If the units of (x,y) are not in meters you will have to adjust $xyfact.
# New lines will be started if the distance between two consecutive points
# will be greater than $maxdist. You have to adjust $maxdist accordingly.
# If you only want one line connecting all points use $maxdist=0.
#
# options:
#
# -maxdist=#
# -xyfact=#
# -grd
#
#-----------------------------------------------------------------

$xyfact = 1. unless $xyfact;    # scale x/y values
$maxdist = 0. unless $maxdist;	# new line if $dist > $maxdist (0 for one line)
$grd = 0. unless $grd;		# input file is in grd format

#-----------------------------------------------------------------

while(<>) {

  chomp;
  next if /^\s*$/;              # empty line
  s/,/ /g;                      # change , to space
  s/;/ /g;                      # change ; to space

  @f = split;

  if( $grd ) {
    $x = $f[3];
    $y = $f[4];
  } else {
    $x = $f[0];
    $y = $f[1];
  }

  $x *= $xyfact;
  $y *= $xyfact;
  $dist = make_dist($x,$y);
  next if $dist == 0.;			# skip non-unique points

  if ( $maxdist > 0 and $dist > $maxdist ) {
    make_line(@nodes);
    @nodes = ();
  }

  make_node($x,$y);
  push(@nodes,$nnode);
}

make_line(@nodes);

#------------------------------------------
# subroutines
#------------------------------------------

sub make_node {

  my ($x,$y) = @_;

  my $ntype = 1;

  $nnode++;

  print "1 $nnode $ntype $x $y\n";
}

sub make_line {

  my @nodes = @_;

  my $ltype = 1;
  my $ntot = @nodes;
  my $i = 0;

  return if $ntot <= 1;

  $nline++;

  print "\n";
  print "3 $nline $ltype $ntot\n";
  foreach (@nodes) {
    print " $_";
    $i++;
    print "\n" if $i%10 == 0;
  }

  print "\n";
}

sub make_dist {

  my ($x,$y) = @_;

  if( not defined $xold ) {	# only for first call
    $xold = $x;
    $yold = $y;
    if( $maxdist > 0 ) {
      return $maxdist/2.;
    } else {
      return 1.;
    }
  }

  my $dx = $x - $xold;
  my $dy = $y - $yold;
  my $dist = $dx*$dx + $dy*$dy;

  $xold = $x;
  $yold = $y;

  return sqrt($dist);
}

#------------------------------------------

