#!/usr/bin/perl -sw
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# usage: add_date.pl [options] -date0=date0 file
#
# adds column to time series where time is in first column
# if no date0 is specified the time column is assumed to be absolute time
# else the time is relative to date0
#
#-------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use strict;
use date;

$::maxcol = 0 unless $::maxcol;
$::date0 = 0 unless $::date0;
$::h = 0 unless $::h;
$::help = 0 unless $::help;

my $date = new date;

my $file = $ARGV[0];
#print STDERR "reading from file $file\n";

if( not defined $file or $::h or $::help ) {
  Usage();
  exit 1;
}

my @res;
my $absolute = 1;
$absolute = 0 if $::date0;
die "absolute time not yet implemented...\n" if $absolute;

$::date0 = "\'$::date0\'";	#must add '' for conversion routine
$date->init_date($::date0);

while(<>) {

  chomp;
  my @f = split;
  my $time = $f[0];

  (@res) = $date->convert_from_it($time);
  (@res) = format_numbers(@res);
  my $dateline = $date->format_time_date(@res);

  #print STDERR "$time  ->  $line\n";

  if( $::maxcol > 0 ) {
    @f = @f[0 .. $::maxcol];	#add extra time column
    $_ = join(" ",@f);
  }
  $_ .= "  $dateline";
  print "$_\n";
}


#  } else {
#    $line = join(" ",@res);

#print "$line\n";

#--------------------------------------

sub Usage {

  print STDERR "Usage: add_date.pl -date0=date0 file\n";
  print STDERR "  adds date column to time series file\n";
  print STDERR "  maxcol    only print up to maxcol columns\n";
  print STDERR "  date0     reference date (YYYY-MM-DD::hh:mm:ss)\n";
  print STDERR "            date may be shortened, e.g. 1997 or 1997-01-01\n";
}

#--------------------------------------

sub format_numbers {

  foreach my $n (@_) {
    if( $n < 1 ) {
      $n = "00";
    } elsif( $n < 10) {
      $n = "0$n";
    }
  }

  return @_;
}

