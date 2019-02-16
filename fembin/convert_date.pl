#!/usr/bin/perl -s
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# usage: convert_date.pl [options] year0 {it|year month day hour min sec}
#
# possible options: -it2date -date2it
#
#-------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use date;

my $date = new date;

$year0 = shift(@ARGV);
@date=@ARGV;

if( not defined $date[0] or $h or $help ) {
  Usage();
  exit 1;
}

if( not $it2date and not $date2it ) {	# no option given
  if( defined $date[1] ) {		# more than one value given
    $date2it = 1;
  } else {
    $it2date = 1;
  }
}

#print "|$it2date|$date2it|\n";

$date->init_year($year0);

if( $it2date ) {
  (@res) = $date->convert_from_it(@date);
  (@res) = format_numbers(@res);
  if( $format ) {
    $line = $date->format_time_date(@res);
  } else {
    $line = join(" ",@res);
  }
} else {
  $date[1] = 1 unless $date[1];		#no month given
  $date[2] = 1 unless $date[2];		#no day given
  $line = $date->convert_to_it(@date);
}

print "$line\n";

#--------------------------------------

sub Usage {

  print STDERR "Usage: convert_date.pl [-it2date|-date2it] year0 {it|date}\n";
  print STDERR "   -it2date:  convert from seconds to date\n";
  print STDERR "   -date2it:  convert from date to seconds\n";
  print STDERR "   -format:   format date\n";
  print STDERR "   year0:     reference year\n";
  print STDERR "   date:      year month day [hour [min [sec]]]\n";
  print STDERR "   it:        time in seconds from 1.1 of year0\n";
  print STDERR " if for {it|date} only one value is given,\n";
  print STDERR " -it2date is the default, else -date2it is used\n";
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

