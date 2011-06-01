#!/usr/bin/perl -s
#
# converts seconds from given year to human readable format
#
#
# Usage: convert_date.pl -year=# ts-file
#
#----------------------------------------------------------

use lib ("$ENV{HOME}/fem/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
use date;

$def_year = 0;
$def_year = $year if( $year );	
unless( $def_year ) {
  die "Must give reference year. Use -year=2005 or similar to set year.\n";
}

my $date = new date;
$date->init_year($def_year);

while(<>) {

  chomp;
  s/^\s+//;
  @f = split;

  $it = shift(@f);
  $rest = join(" ",@f);

  ($year,$month,$day,$hour,$min,$sec) = $date->convert_from_it($it);
  ($month,$day,$hour,$min,$sec) = format_numbers($month,$day,$hour,$min,$sec);

  #print "$year/$month/$day $hour:$min:$sec  $rest\n";
  print "$day/$month/$year $hour:$min:$sec  $rest\n";
}

print STDERR "Using $def_year as reference year\n";

#------------------------------------------------------------

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
