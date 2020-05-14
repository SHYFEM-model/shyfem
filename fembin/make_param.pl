#!/usr/bin/perl -w

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use strict;

#---------------------------------------------- initialization

my $debug = 0;

my $dir = shift;
my $file = shift;

my $parfile = "$dir/$file";
my @newlines = ();
my $ichange = 0;
my %dims = ();

#---------------------------------------------- read environment

foreach my $param (sort keys %ENV) {
  my $val = $ENV{$param};
  if( $param =~ /^\w+DIM$/ ) {
    my $name = lc($param);
    $dims{$name} = $val;
    print "$param ($name) = $val\n" if $debug;
  }
}

#---------------------------------------------- read param file

open(FILE,"$parfile") || die "Cannot open file: $parfile\n";
my @lines = <FILE>;
close(FILE);

#---------------------------------------------- check differences

foreach my $line (@lines) {
  if( $line =~ /^\s+parameter\s*\(\s*(\w+dim)\s*=\s*(\d+)\s*\)\s*(!.*)?$/i ) {
    my $param = $1;
    my $val = $2;
    print "parameter found: $param = $val\n" if $debug;
    if( $dims{$param} ) {
      my $newval = $dims{$param};
      print "parameter is given: $newval\n" if $debug;
      if( $newval != $val ) {
	$ichange++;
	$line = "\tparameter ( $param = $newval )\n";
      }
    }
  }
  push(@newlines,$line);
}

#---------------------------------------------- if different write param file

exit 0 unless $ichange;

print "parameters have changed ($ichange changes): writing $parfile\n";

rename($parfile,"$parfile.bak");

open(FILE,">$parfile") || die "Cannot open file: $parfile\n";
foreach my $line (@newlines) {
  print FILE "$line";
}
close(FILE);

#---------------------------------------------- end of routine

