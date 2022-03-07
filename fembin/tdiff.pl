#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# computes time difference
#
# Usage: timediff time1 time2
#
# where time is given in standard date format: YYYY-MM-DD[::hh:mm:ss]
#
#------------------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use date;
 
my $date = new date;

my $time1 = $ARGV[0];
my $time2 = $ARGV[1];

unless( $time2 ) {
  print STDERR "Usage: tdiff.pl time1 time2\n";
  print STDERR "  computes difference between times (time2-time1) in seconds\n";
  print STDERR "  time must be in format YYYY-MM-DD[::hh:mm:ss]\n";
  exit 0
}

print "time1: $time1\n";
print "time2: $time2\n";

my $atime1 = $date->unformat_abs($time1);
my $atime2 = $date->unformat_abs($time2);

my $difftime = $atime2 - $atime1;

print "diff: $difftime\n";

