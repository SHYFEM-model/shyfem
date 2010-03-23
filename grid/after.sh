#!/usr/bin/perl -s

# extracts entry from revision history on a specified date
#	and formats it to be included in gridlog
# can give -date on command line
# example:      "after -date 20-Oct-97 *.c"

%months = (
		  "..." => "0" 
		, "Jan" => "1"
		, "Feb" => "2"
		, "Mar" => "3"
		, "Apr" => "4"
		, "May" => "5"
		, "Jun" => "6"
		, "Jul" => "7"
		, "Aug" => "8"
		, "Sep" => "9"
		, "Oct" => "10"
		, "Nov" => "11"
		, "Dec" => "12"
	  );

if( $date ) {
  $date = shift;
  $afterdate = &getdate($date);
} else {
  &usage;
}

while(<>) {
  $revhist = 0 if( /^ \*\s+\*\s*\n$/ );
  &elabline if( $revhist );
  $revhist = 1 if( /^ \* Revision History:\s+\*\s*\n$/ );
}


############ subroutines ##############

sub elabline {

  if( /^ \* (.+-.+-.+): (.*)\s*\*\s*\n$/ ) {	#date line
	$date = $1;
	$info = $2;
	$cdate = &getdate($date);
  } elsif( /^ \*\s+(.+)\s*\*\s*\n$/ ) {		#conti line
	$info = $1;
  } else {					#error
	die "Format error: $_";
  }

  if( $cdate >= $afterdate ) {
	&printline;
  }
}
	
sub getdate {

  my $date = $_[0];
  my ($day,$month,$year,$new);

  if( $date =~ /^(.+)-(.+)-(.+)$/ ) {
    $day = $1;
    $month = $2;
    $year = $3;

    $day = 0 if( $day eq '..' );
    $year += 1900 if( $year < 100 );

    $month = $months{$month};
    if( $month eq "" ) {
	print STDERR "Error in month: $date\n";
	$month = 0;
    }
    if( $day < 0 || 31 < $day ) {
	print STDERR "Error in day: $date\n";
	$day = 0;
    }
    if( $year < 1970 || 2100 < $year ) {
	print STDERR "Error in year: $date\n";
	$year = 0;
    }

    $new = $day + 100 * $month + 10000 * $year;
  } else {
    print STDERR "Error in date: $date\n";
    $new = 0;
  }

  return $new;
}

sub printline {
  s/^ \*//;
  s/\s*\*\s*\n$//;
  if( $ARGV ne $oldfile ) {
    print "$ARGV:\n";
    $oldfile = $ARGV;
  }
  print "$_\n";
}

sub usage {
  print "Usage: after -date \"date-string\" file(s)\n";
  print "   (example: after -date 20-Oct-97 *.c)\n";
  exit;
}
