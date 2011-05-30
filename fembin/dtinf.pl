#!/usr/bin/perl -w
#
# analyses INF file and writes info about internal time steps
#
#-------------------------------------------------------------

use strict;

my $itold = -1;
my $out0;
my $acum = 0;
my $ripet = 0;
my $all = 0;
my $all0 = 0;
my $elem0 = 0;

while(<>) {

  chomp;

  if( /^newlnk:/ ) {
    my @f = split;
    my $it = $f[1];
    my $out = $f[11];
    my $elem = $f[3];

    #print "$it $out    $itold\n";

    if( $it != $itold ) {
      if( $itold != -1 ) {
        print "$itold   $ripet  $acum $elem0\n";	#ripet and number of elements
      }
      $all += $ripet + 1;
      $all0++;
      $acum = 0;
      $ripet = 0;
      $out0 = $out;
    } else {
      $ripet++;
      $acum = $out - $out0;
    }
    $itold = $it;
    $elem0 = $elem;
  }
}

print "$itold   $ripet  $acum  $elem0\n";

$all += $ripet + 1;
$all0++;

my $tot = $all/$all0;

# $ripet	number of repetitions of time step (minimum 0)
# $acum		number of elements excluded during iterations
# $all		number of internal time steps (minimum 1)
# $all0		number of external time steps
# $tot		ratio of time steps (best minimum 1)

print STDERR "$all  $all0  $tot\n";

