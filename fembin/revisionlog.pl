#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/fembin","$ENV{HOME}/shyfem/fembin");

require "revision_getdate.pl";

# gets revision log from header of fortran or c file

%months = (
		 "Jan"	=>	1
		,"Feb"	=>	2
		,"Mar"	=>	3
		,"Apr"	=>	4
		,"May"	=>	5
		,"Jun"	=>	6
		,"Jul"	=>	7
		,"Aug"	=>	8
		,"Sep"	=>	9
		,"Oct"	=>	10
		,"Nov"	=>	11
		,"Dec"	=>	12
	  );

$debug = 0;

$after = 0;
$befor = 30000000;
$noname = 0;
$sepname = "";
$check = 0;

for(;;) {
  $option = $ARGV[0];
  if( $option eq "-check" ) { shift; $check = 1; next; }
  if( $option eq "-after" ) { shift; $after = shift; next; }
  if( $option eq "-befor" ) { shift; $befor = shift; next; }
  if( $option eq "-file" ) { shift; $file = shift; next; }
  if( $option eq "-noname" ) { shift; $noname = 1; next; }
  if( $option eq "-sepname" ) { shift; $sepname = "\n"; next; }
  last;
}

$savefile = $file;

if( $check ) {
  #print STDERR "checking $file\n";
  revisionlog_fortran_check();
  exit 0;
}

while(<>) {

  if( /^[cC\!]\s+revision log\s+:\s*\n$/ ) {
	revisionlog_fortran();
	last;
  } elsif( /^\s+\*\s+Revision History:\s+\*\s*\n$/ ) {
	revisionlog_c();
	last;
  }

}

sub revisionlog_c {

  my $olddate = $befor;		#in c files order is inverse

  while(<>) {
    if( /^\s+\*\s+\*\s*\n$/ ) {			#only white space
	last;
    }
    if( /^\s+\\\*+\/\s*\n$/ ) {			# \****...***/
	last;
    }
    $date = &getdate($date);
    print STDERR "     $after-$date-$befor\n$_" if $debug;
    if( $date && $after <= $date && $date <= $befor ) {
      if( $file && not $noname ) {
	print "${sepname} * $file:\n${sepname}";
	$file = "";
      }
      print;
    }
    if( $date > $olddate ) {
	print STDERR "*** error in date: $file  $olddate $date\n";
	print STDERR "       $_";
    }
    $olddate = $date;
  }
}

sub revisionlog_fortran {

  my $olddate = 0;

  while(<>) {
    if( /^[cC\!]\s+(.*)\s+:\s*\n$/ ) {		#new keyword
	last;
    }
    if( /^\s*\n$/ ) {				#only white space
	last;
    }
    if( /^[cC\!]\*+\s*\n$/ ) {			# c****
	last;
    }
    if( /^[cC\!]-+\s*\n$/ ) {			# !----
	last;
    }
    $date = &getdate($date);
    print STDERR "     $after-$date-$befor $_" if $debug;
    if( $date && $after <= $date && $date <= $befor ) {
      if( $file && not $noname ) {
	print "${sepname}c $file:\n${sepname}";
	$file = "";
      }
      print;
    }
    if( $date < $olddate ) {
	print STDERR "*** error in date: $file  $olddate $date\n";
	print STDERR "       $_";
    }
    $olddate = $date;
  }
}

sub revisionlog_fortran_check {

  my $rev_log_found = 0;
  my $ndate = 0;
  my $debug = 0;
  my $olddate = 0;

  while(<>) {

    if( /^\s*\n$/ ) {				#only white space
	next;
    }
    if( /^\s+\S$/ ) {				# first statement
	last;
    }
    if( /^[cC\!]\s+revision log\s+:\s*\n$/ ) {
	$rev_log_found = 1;
	next;
    }

    $date = &getdate($date);

    if( $date < $olddate ) {
	print STDERR "*** error in date: $file  $olddate $date\n";
	print STDERR "       $_";
    }
    $olddate = $date;

    if( $date ) {
      $ndate++;
      print if $debug;
    }
  }

  print "$file:  $rev_log_found  $ndate";
  print "   ********" if $rev_log_found == 0 and $ndate > 0;
  print "\n"
}

#---------------------------------------------------------------

sub getdate {

  my $olddate = shift;

  my $date = get_date($_,$savefile);

  $date = $olddate unless $date;

  return $date;
}

#---------------------------------------------------------------

