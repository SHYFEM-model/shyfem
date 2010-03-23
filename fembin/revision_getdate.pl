#!/usr/bin/perl

# gets date of revision log from header of fortran or c file

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

sub get_date {

  my $line = shift;
  my $date;

  if( $line =~ /^[cC\!]\s+(\d+)\.(\d+)\.(\d+)\s+/ ) {

	my $day = $1;
	my $mon = $2;
	my $y = $3;

	die "*** Error in date: ($savefile)\n$_" if( $day < 0 || $day > 31 );
	die "*** Error in date: ($savefile)\n$_" if( $mon < 0 || $mon > 12 );
	die "*** Error in date: ($savefile)\n$_" if( $y < 1000 || $y > 2500 );

	$date = 10000 * $y + 100 * $mon + $day;

  } elsif( $line =~ /^\s+\*\s+(\d+)-(\w+)-(\d+):\s+/ ) {

	my $day = $1;
	my $mon = $2;
	my $y = $3;

	if( $mon =~ /^[a-zA-Z]+$/ ) {
	  my $month = $months{$mon};
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

	$date = 0;

  }

  return $date;
}

sub is_date {

  my $line = shift;

  return get_date($line);
}

sub is_file_name {

  my $line = shift;

  if( $line =~ /^[cC\!]\s+(\S+):\s*\n$/ ) {
    return $1;
  } elsif( $line =~ /^\s*\*\s+(\S+):\s*\*?\s*\n$/ ) {
    return $1;
  } else {
    return "";
  }
}

sub is_empty_comment {

  my $line = shift;

  if( $line =~ /^[cC\!]\s*\n$/ ) {
    return 1;
  } elsif( $line =~ /^\s*\*\s+\*\s*\n$/ ) {
    return 1;
  } else {
    return 0;
  }
}

1;

