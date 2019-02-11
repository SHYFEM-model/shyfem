#!/usr/bin/perl -ws
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
#-----------------------------------

use strict;

$::ext= "" unless $::ext;
$::quiet= "" unless $::quiet;

my $file = $ARGV[0];

unless( $file ) {
  die "Usage: check_copyright.pl [-ext extension] file\n";
}

handle_file();

print "no copyright found\n" unless $::quiet;
exit 1;

#-----------------------------------

sub handle_file
{
  my $debug = 0;	# set to > 0 for debug
  my $shyfem = 0;
  my $copy = 0;
  my $contrib = 0;

  while(<>) {
    if( /Copyright \(C\) (.*)/ ) {
      $copy = $1;
    }
    if( /This file is part of SHYFEM\./ ) {
      $shyfem = 1;
    }
    if( $copy and $shyfem ) {	# we have found copyright
      print "$copy\n" unless $::quiet;
      exit 0;
    }
  }
}

#-----------------------------------

