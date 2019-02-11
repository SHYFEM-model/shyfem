#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

sub days {
  my $day = shift;
  $day = lc($day);
  return $days{$day};
}

%days = (
		 "mon"	=>	"1"
		,"tue"	=>	"2"
		,"wed"	=>	"3"
		,"thu"	=>	"4"
		,"fri"	=>	"5"
		,"sat"	=>	"6"
		,"sun"	=>	"7"
	);

sub months {
  my $month = shift;
  $month = lc($month);
  return $months{$month};
}

%months = (
		 "jan"	=>	"01"
		,"feb"	=>	"02"
		,"mar"	=>	"03"
		,"apr"	=>	"04"
		,"may"	=>	"05"
		,"jun"	=>	"06"
		,"jul"	=>	"07"
		,"aug"	=>	"08"
		,"sep"	=>	"09"
		,"oct"	=>	"10"
		,"nov"	=>	"11"
		,"dec"	=>	"12"
	  );


sub make_date {

  my $date = shift;

  my $save = $date;
  my $unknown = "00000000";

  my ($day,$month,$year);

  return $unknown unless $date;

  my @f = split(/\s+/,$date);
  my $nf = @f;

  if( $f[0] =~ /([a-zA-Z]+),/ ) { 	# Fri, 19 Nov 1999 14:21:21 +0100
    $day = $1;
    if( &days($day) ) {
      $day = $f[1];
      $month = &months($f[2]);
      $year = $f[3];
    } else {
      print STDERR "Unknown format of date (1): $date\n";
      return $unknown;
    }
  } elsif( $f[0] =~ /([a-zA-Z]+)/ ) { # Fri Jan 20 10:17:41 1995
    $day = $1;				# Thu Dec 23 17:30:10 MET 1999
    if( &days($day) ) {
      $day = $f[2];
      $month = &months($f[1]);
      $year = $f[$nf-1];
    } else {
      print STDERR "Unknown format of date (2): $date\n";
      return $unknown;
    }
  } elsif( $f[1] =~ /([a-zA-Z]+)/ ) { 	# 12 Apr 1999 17:42:24 +0900
    $month = $1;
    if( &months($month) ) {
      $day = $f[0];
      $month = &months($f[1]);
      $year = $f[2];
    } else {
      print STDERR "Unknown format of date (3): $date\n";
      return $unknown;
    }
  } else {
      print STDERR "Unknown format of date (9): $date\n";
      return $unknown;
  }

  if( $day < 1 || $day > 31 ) {
      print STDERR "Error in date: $date  (1)  ($day/$month/$year)\n";
      return $unknown;
  }
  if( $month < 1 || $month > 12 ) {
      print STDERR "Error in date: $date  (2)  ($day/$month/$year)\n";
      return $unknown;
  }
  if( $year < 1 || $year > 100 && $year < 1980 || $year > 2030 ) {
      print STDERR "Error in date: $date  (3)  ($day/$month/$year)\n";
      return $unknown;
  }

  $day += 0;		#convert to numeric ("05" -> 5)
  $day = "0" . $day if $day < 10;
  $year += 1900 if $year < 1000;

  #$date = "$f[1]-$f[2]-$f[3]";
  $date = "$year$month$day";

  if( $date < 19800000 || $date > 20300000 ) {
      print STDERR "Error in date conversion: $save  ->  $date\n";
      return $unknown;
  }

  return $date;
}

sub format_date
{
  my $date = shift;

  $aux = int(0.5 + $date/100);
  $day = $date - 100*$aux;
  $date = $aux;

  $aux = int(0.5 + $date/100);
  $month = $date - 100*$aux;
  $date = $aux;

  $year = $date;

  return "$day-$month-$year";
}

#########
1;
#########


