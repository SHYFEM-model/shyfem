#!/usr/bin/perl -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# converts time column in seconds given a format
#
#--------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
use date;
use strict;

#-------------------------------------------------
#-------------------------------------------------
#-------------------------------------------------
# please set time/date format and reference year
#-------------------------------------------------
#-------------------------------------------------
#-------------------------------------------------

# examples for format:
# DD-MM-YYYY hh:mm:ss
# MM/DD/YYYY hh:mm:ss

my $format = "DD/MM/YYYY hh:mm";	# format for time/date column
my $year0 = 2014;			# reference year from when to count
my $add_formatted_date = 1;		# append date to data columns

#-------------------------------------------------
#-------------------------------------------------
#-------------------------------------------------
# do not change anything beyond here
#-------------------------------------------------
#-------------------------------------------------
#-------------------------------------------------

$::debug = 0;
$::add_formatted_date = $add_formatted_date;

%::descrp = 
	( "D" => '$day'
	, "M" => '$month'
	, "Y" => '$year'
	, "h" => '$hour'
	, "m" => '$min'
	, "s" => '$sec'
	);

my $date = new date;

$date->init_year($year0);             #set reference year

my ($pattern,$eline) = assemble_pattern($format);

while(<>) {

  chomp;
  next if /^\s*$/;

  my @matches = (m/^\s+$pattern\s+(.+)$/g);
  my $n = @matches;
  if( $n ) {
    my $rest = pop(@matches);
    convert($date,\@matches,$eline,$rest);
  } else {
    die "Cannot parse time: $_\n";
  }

}

#---------------------------------------------

sub convert
{
  my ($date,$fields,$eline,$rest) = @_;

  my ($year,$month,$day,$hour,$min,$sec) = (0,0,0,0,0,0);

  my $line = join(",",@$fields);
  my $final = "$eline = ($line);";

  eval "$eline = ($line);";

  if( $::debug ) {
    print STDERR "$line ... $rest\n";
    print STDERR "eval: $final\n";
    print STDERR "$year $month $day $hour $min $sec\n";
  }

  my $it = $date->convert_to_it($year,$month,$day,$hour,$min,$sec);
  my $formatted = $date->format_time_date($year,$month,$day,$hour,$min,$sec);

  if( $::add_formatted_date ) {
    print "$it  $rest  $formatted\n";
  } else {
    print "$it  $rest\n";
  }
}

sub assemble_pattern
{
  my $format = shift;

  my $pattern = "";
  my $what = [];
  my @array = split(//,$format);

  my $last = "";

  foreach my $c (@array) {

    if( $c =~ /(\w+)/ ) {		#letter ... stands for digit
      my $char = $1;
      next if $char eq $last;
      die "error in format: $format\n" if $last;
      push(@$what,$char);
      $pattern .= "(\\d+)";
      $last = $char;
    } elsif( $c =~ /\s+/ ) {		#white space
      $pattern .= "\\s+";
      $last = "";
    } else {
      $pattern .= "\\$c";
      $last = "";
    }
  }

  my @eline = ();
  foreach my $item (@$what) {
    my $trans = $::descrp{$item};
    unless( $trans ) {
      die "do not understand format: $item\n"
    }
    push(@eline,$trans);
  }
  my $eline = "(" . join(",",@eline) . ")";

  my $line = join(" ",@$what);
  if( $::debug ) {
    print STDERR "pattern: $pattern\n";
    print STDERR "what: $line\n";
    print STDERR "eline: $eline\n";
  }

  return ($pattern,$eline);
}

