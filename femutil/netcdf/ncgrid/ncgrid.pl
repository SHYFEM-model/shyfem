#!/usr/bin/perl

#--------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#--------------------------------------------------------------------------
#
# parses nc header for information on lat/lon
#
#---------------------------------------------------

$var = shift;
print STDERR "parsing variable $var\n";

while(<>) {

  chomp;

  if( /^dimensions:/ ) {
    $in_dimensions = 1;
    next;
  }
  if( /^variables:/ ) {
    $in_dimensions = 0;
    $in_variables = 1;
    next;
  }
  if( /global attributes:/ ) {
    print STDERR "... in global\n";
    $in_variables = 0;
    next;
  }
  if( /$var =/ ) {
    print STDERR "... in data\n";
    my $xsize = $::dims{x};
    my $ysize = $::dims{y};
    print STDERR "... size: $xsize $ysize\n";
    print " $xsize $ysize\n";
    $in_data = 1;
    next;
  }

  if( $in_dimensions ) {
    s/\s*;.*$//;
    s/^\s+//;
    insert_dimension($_);
    print STDERR "dimension: $_\n";
  }
  if( $in_variables ) {
    s/\s*;\s*$//;
    if( /double $var(.*)$/ ) {
      print STDERR "variable $var: $1\n";
    }
  }
  if( $in_data ) {
    s/^\}$//;
    s/\s*;\s*$//;
    print "$_\n";
  }
}

sub insert_dimension
{
  my $dim = shift;

  my @f = split(/\s*=\s*/,$dim);

  my $what = $f[0];
  my $size = $f[1];
  $::dims{$what} = $size;
}

