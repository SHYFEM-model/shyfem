#!/usr/bin/perl -ws
#
#-----------------------------------

use strict;

$::ext= "" unless $::ext;

my $file = $ARGV[0];

unless( $file ) {
  die "Usage: check_copyright.pl [-ext extension] file\n";
}

handle_file();

print "no copyright found\n";
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
      print "$copy\n";
      exit 0;
    }
  }
}

#-----------------------------------

