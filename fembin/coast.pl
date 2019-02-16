#!/usr/bin/perl -s
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
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
# -maxdist=#		maximum distance for connecting in line
# -xyfact=#		apply factor to x/y values
# -grd			input file is in grd format
# -invert		invert x/y (e.g., if given in lat/lon)
# -close		close created lines
#
#-----------------------------------------------------------------

$maxdist = 0. unless $maxdist;	# new line if $dist > $maxdist (0 for one line)
$xyfact = 1. unless $xyfact;    # scale x/y values
$grd = 0. unless $grd;		# input file is in grd format
$invert = 0. unless $invert;	# invert x/y
$close = 0. unless $close;	# close created lines
$help = 0. unless $help;	# write help
$h = 0. unless $h;	# write help

if( $help or $h ) {
  print STDERR "coast.pl [-h|-help] [-options] file(s)\n";
  print STDERR "   Takes single points and creates lines.\n";
  print STDERR "   options:\n";
  print STDERR "     -h|-help     this help screen\n";
  print STDERR "     -grd         input file is in grd format\n";
  print STDERR "     -invert      invert x/y\n";
  print STDERR "     -close       close created lines\n";
  print STDERR "     -maxdist=#   new line if $dist > # (0 for one line)\n";
  print STDERR "     -xyfact=#    apply factor to x/y values\n";
  exit 1;
}

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

  ($x,$y) = invert($x,$y) if $invert;

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
  $x[$nnode] = $x;
  $y[$nnode] = $y;

  #print "1 $nnode $ntype $x $y\n";
}

sub make_line {

  my @nodes = @_;

  my $ntype = 1;
  my $ltype = 1;
  my $ntot = @nodes;
  my $i = 0;

  return if $ntot <= 1;

  $nline++;

  my $closed = 0;
  if( equal_node($nodes[0],$nodes[$ntot-1]) ) {
    pop(@nodes);
    push(@nodes,$nodes[0]);
    $closed = 1;
  } elsif( $close ) {
    push(@nodes,$nodes[0]);
    $ntot++;
    $closed = 1;
  }

  my $nlast = pop(@nodes) if $closed;
  print "\n";
  foreach my $n (@nodes) {
    my $x = $x[$n];
    my $y = $y[$n];
    print "1 $n $ntype $x $y\n";
  }
  push(@nodes,$nlast) if $closed;

  print "\n";
  print "3 $nline $ltype $ntot\n";
  foreach (@nodes) {
    print " $_";
    $i++;
    print "\n" if $i%10 == 0;
  }

  print "\n";
}

sub equal_node {

  my ($i1,$i2) = @_;

  if( $x[$i1] == $x[$i2] and $y[$i1] == $y[$i2] ) {
    return 1;
  } else {
    return 0;
  }
}

sub invert {
  return ($_[1],$_[0]);
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

