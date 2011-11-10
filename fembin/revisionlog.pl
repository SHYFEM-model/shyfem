#!/usr/bin/perl

$FEMDIR=$ENV{SHYFEMDIR}?$ENV{SHYFEMDIR}:$ENV{HOME}/shyfem;

push(@INC,"$FEMDIR/fembin");

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
  }
}

sub revisionlog_fortran {

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
    $date = &getdate($date);
    print STDERR "     $after-$date-$befor $_" if $debug;
    if( $date && $after <= $date && $date <= $befor ) {
      if( $file && not $noname ) {
	print "${sepname}c $file:\n${sepname}";
	$file = "";
      }
      print;
    }
  }
}

sub revisionlog_fortran_check {

  my $rev_log_found = 0;
  my $ndate = 0;
  my $debug = 0;

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
}


#---------------------------------------------------------------

sub getdate {

  my $olddate = shift;

  my $date = get_date($_);

  $date = $olddate unless $date;

  return $date;
}

sub getdate0 {

  my $olddate = shift;
  my $date;

  if( /^[cC\!]\s+(\d+)\.(\d+)\.(\d+)\s+/ ) {

	my $day = $1;
	my $mon = $2;
	my $y = $3;

	die "*** Error in date: ($savefile)\n$_" if( $day < 0 || $day > 31 );
	die "*** Error in date: ($savefile)\n$_" if( $mon < 0 || $mon > 12 );
	die "*** Error in date: ($savefile)\n$_" if( $y < 1000 || $y > 2500 );

	$date = 10000 * $y + 100 * $mon + $day;

  } elsif( /^\s+\*\s+(\d+)-(\w+)-(\d+):\s+/ ) {

	my $day = $1;
	my $mon = $2;
	my $y = $3;

	if( $mon =~ /^[a-zA-Z]+$/ ) {
	  $month = $months{$mon};
	  die "*** Error in month: ($savefile)\n$_" unless $month;
	  $mon = $month;
	}

	if( $y < 100 && $y > 70 ) {
	  $y += 1900;
	}

	die "*** Error in date: ($savefile)\n$_" if( $day < 0 || $day > 31 );
	die "*** Error in date: ($savefile)\n$_" if( $mon < 0 || $mon > 12 );
	die "*** Error in date: ($savefile)\n$_" if( $y < 1000 || $y > 2500 );

	$date = 10000 * $y + 100 * $mon + $day;

  } else { 			# use old date

	$date = $olddate;

  }

  return $date;
}
