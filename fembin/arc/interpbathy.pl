#!/usr/bin/perl -w -s
#
# interpolates bathymetry
#
# substitutes bastreat

use lib "$ENV{HOME}/lib/perl";
use grd;
use strict;

# options --------------------------------------------------------------

$::help = 1 if $::help;
$::h = 1 if $::h;
$::all = 0 unless $::all;	#if set interpolate also in elements with depth
$::fact = 1 if $::fact;

if( $::h or $::help ) {
  FullUsage();
  exit 0;
} elsif( not $ARGV[1] ) {
  Usage();
  exit 1;
}

# main -----------------------------------------------------------------

my $grid = new grd;
my $bathy = new grd;
my $gfile = $ARGV[0];	#grid file
my $bfile = $ARGV[1];	#bathymetry file

$grid->readgrd($gfile);
$bathy->readgrd($bfile);

print STDERR "starting interpolation...\n";

delete_depth($grid) if $::all;
make_interp($grid,$bathy,$::fact);

$grid->writegrd("newbathy.grd");;

#-----------------------------------------------------------------------

sub FullUsage {
  Usage();
  print STDERR "  Interpolates bathymetry from bathy into grid.\n";
  print STDERR "  \n";
  print STDERR "  grid must be a finite element grid\n";
  print STDERR "  bathy must contain single nodes with depth values\n";
  print STDERR "  \n";
  print STDERR "  Only elements without depth values are interpolated.\n";
  print STDERR "  The routine uses the element area for the standard\n";
  print STDERR "    deviation (sigma) of the gauss curve used in\n";
  print STDERR "    interpolation of the depth values. The value of fact\n";
  print STDERR "    is used to set a different size of sigma which is\n";
  print STDERR "    taken as sogma_0 = fact*area.\n";
  print STDERR "  If after an iteration elements are still without depth,\n";
  print STDERR "    the procedure is repeated with fact=2*fact.\n";
  print STDERR "  \n";
  print STDERR "  options:\n";
  print STDERR "  -h|-help     this help screen\n";
  print STDERR "  -all         interpolate also in elements with depth\n";
  print STDERR "  -fact=fact   use fact as start value (default=1)\n";
}

sub Usage {
  print STDERR "Usage: interpbathy.pl [-h|-help] [-options] grid bathy\n";
}

# do interpolation ---------------------------------------------------------

sub make_interp
{
  my $grid = shift;
  my $bathy = shift;
  my $fact = shift;

  $fact = 1. if $fact <= 0.;

  my $ntot = 0.;
  my @elems = $grid->get_elem_list();
  my @nodes = $bathy->get_node_list();
  my $nel = @elems;

  while( $ntot < $nel ) {
    print STDERR "iteration $fact\n";
    $ntot = 0.;
    my $nnew = 0.;
    foreach my $elem (@elems) {
      my $depth = $elem->{h};

      if( defined $depth ) {
        $ntot++;
      } else {
	$nnew += interpol($grid,$bathy,$elem,$fact,\@nodes);
      }
      print STDERR "$fact $ntot $nnew\n";
    }
    my $ntotal = $ntot + $nnew;
    my $nmiss = $nel - $ntotal;
    print STDERR "$fact  $nmiss  $ntot  $ntotal  $nel\n";
    $fact *= 2.;
  }

}

# interpolate in one element ------------------------------------------------

sub interpol
{
  my ($grid,$bathy,$elem,$fact,$nodes) = @_;

  my $number = $elem->{number};
  my ($area,$xt,$yt) = area_ct($grid,$elem);
  my $r2max = $fact * $area;		#maximum radius
  my $sigma2 = $fact * $area;		#standard deviation also grows

  my $depth = 0;
  my $weight = 0;
  my $n = 0;

  foreach my $node (@$nodes) {
    my $x = $node->{x};
    my $y = $node->{y};
    my $d = $node->{h};
    my $r2 = ($x-$xt)*($x-$xt) + ($y-$yt)*($y-$yt);
    if( $r2 <= $r2max ) {
      my $w = exp(-$r2/$sigma2);
      $depth += $d * $w;
      $weight += $w;
      $n++;
    }
  }

  if( $n ) {
    if( $weight <= 0. ) {
      die "zero weight from points: $n $weight $r2max $number\n";
    }
    $depth /= $weight;
    $elem->{h} = $depth;
    return 1;
  } else {
    return 0;
  }
}

# deletes depth in all elements ---------------------------------------------

sub delete_depth
{
  my $grid = shift;

  my @elems = $grid->get_elem_list();

  foreach my $elem (@elems) {
    my $depth = $elem->{h};
    if( defined $depth ) {
      undef($elem->{h});
    }
  }
}

# computes area and center point of element (triangle) ----------------------

sub area_ct 
{
  my ($grid,$elem) = @_;

  my @x = ();
  my @y = ();
  my $xt = 0;
  my $yt = 0;

  my $vert = $elem->{vert};

  foreach my $number (@$vert) {
    my $node = $grid->get_node($number);
    my $x = $node->{x};
    my $y = $node->{y};
    push(@x,$x);
    push(@y,$y);
    $xt += $x;
    $yt += $y;
  }

  my $n = @$vert;
  $xt /= $n;
  $yt /= $n;

  if( $n != 3 ) {
    die "area_ct: Cannot compute area for $n > 3\n";
  }
  #for(my $i=0; $i<$n; $i++) {
  #}

  my $area = 0.5* ( ($x[1]-$x[0])*($y[2]-$y[0]) - ($x[2]-$x[0])*($y[1]-$y[0]) );

  return ($area,$xt,$yt);
}

# ---------------------------------------------------------------------------

