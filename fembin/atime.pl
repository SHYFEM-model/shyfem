#!/usr/bin/perl -ws
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# returns absolute time from formatted date/time
#
# date must be given in following format: YYYY-MM-DD::hh:mm:ss
# smaller units can be left out and will be interpolated
# but the first dash after year must exists to distinguish it
# from a normal integer
# so 1999 is not allowed, but 1999- is, as well as 1999-01 or 1999-01-01
#
#-----------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use date;
use strict;

my $debug = 0;

my $line = $ARGV[0];
$line="" unless $line;

my ($date,$time);

my @f = split(/::/,$line);
my $n = @f;

if( $n > 2 ) {
  die "Cannot parse date: $line\n";
} elsif( $n == 2 ) {
  $date = $f[0];
  $time = $f[1];
} elsif( $n == 1 ) {
  $date = $f[0];
  $time = 0;
} elsif( $n == 0 ) {
  print "0\n";
  exit 0;
} 

unless( $date =~ /^\d+-/ ) {	#only number, no date
  print "0\n";
  exit 0;
}

my $pdate = new date;

print STDERR "atime.pl debug: $date  $time\n" if $debug;

my ($hour,$min,$sec) = split(/:/,$time);
my ($year,$mon,$day) = split(/-/,$date);

my $atime = $pdate->convert_to_abs($year,$mon,$day,$hour,$min,$sec);
my $dline = $pdate->format_abs($atime);

print STDERR "atime.pl debug: date parsed as: $dline\n" if $debug;

print "$atime\n";

#----------------------------------------------------

