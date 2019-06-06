#!/usr/bin/perl -ws

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use strict;

$::inverse = 0 unless $::inverse;

$::n = 0;
$::l = 0;
@::nums = ();

my $file = $ARGV[0];

if( $::inverse ) {
  bnd2grd($file);
} else {
  grd2bnd($file);
}

#----------------------------------------------------------

sub grd2bnd
{
  my $file = shift;

  my $grid = new grd;
  $grid->readgrd($file);

  my $lines = $grid->get_lines();

  foreach my $line (values %$lines) {
    my $flag = 1;
    my $vert = $line->{vert};
    foreach my $nnode (@$vert) {
      my $node = $grid->get_node($nnode);
      my $x = $node->{x};
      my $y = $node->{y};
      print "$x $y   $flag\n";
      $flag = 0;
    }
  }
}

#----------------------------------------------------------

sub bnd2grd
{
  my $file = shift;

  while(<>) {
    chomp;
    my ($x,$y,$flag) = split;
    if( $flag ) {
	write_line();
    }
    $::n++;
    print "1 $::n 0 $x $y\n";
    push(@::nums,$::n)
  }
  write_line();
}

#----------------------------------------------------------

sub write_line
{
  my $n = @::nums;

  return if $n == 0;

  $::l++;
  print "3 $::l 0 $n\n";

  my $c = 0;
  foreach my $nn (@::nums) {
    print " $nn";
    $c++;
    print "\n" if( $c%10 == 0 );
  }
  print "\n" if( $c%10 != 0 );

  @::nums = ();
}

#----------------------------------------------------------

