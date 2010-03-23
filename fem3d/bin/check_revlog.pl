#!/usr/bin/perl

$file = $ARGV[0];

while(<>) {
  if ( /^[cC\!]\s+(revision log)\s+:\s*\n$/ ) {
    $rev_found = 1;
    last;
  }
}

if( $rev_found ) {
  @revlog = get_revisionlog_fortran();
} else {
  #print STDERR "missing revision log\n";
  @revlog = ();
}

#print_array(@revlog);

check_revisionlog(@revlog);
#print_dates();

#---------------------------------------------------------------

sub get_revisionlog_fortran {

  my $revlog = ();

  while(<>) {
    last if( /^[cC\!]\s+(.*)\s+:\s*\n$/ ) ;		#new keyword
    last if( /^\s*\n$/ ) ;				#only white space
    last if( /^[cC\!]\*+\s*\n$/ ) ;			#c****
    next if( /^[cC\!]\s*\n$/ ) ;			#c

    chomp;
    push(@revlog,$_);
  }

  return @revlog;
}

sub getdate {

  my $olddate = shift;

  my $date = 0;

  if( $olddate =~ /^(\d+)\.(\d+)\.(\d+)$/ ) {
    my $day = $1;
    my $month = $2;
    my $year = $3;

    $year += 1900 if( $year < 100 );
    $date = 10000*$year + 100*$month + $day;
  } else {
    #die "*** Cannot parse date: $olddate\n";
    print STDERR "*** Cannot parse date: $olddate\n";
  }

  return $date;
}

sub print_array {

  foreach (@_) {
    print "$_\n";
  }
}

sub print_dates {

  my @keys = sort keys %dates;
  foreach my $key (@keys) {
    my $count = $dates{$key};
    print "$key :  $count\n";
  }
}

sub check_revisionlog {

  my $olddate = 0;

  foreach my $line (@_) {
    ($date,$who,$comment) = check_revline($line);
    if( $date ) {
      $date = getdate($date);
      $dates{$date}++;
      #print "$date: $comment\n";
      if( $date < $olddate ) {
	print STDERR "*** Error in date: $olddate $date\n";
      }
      $olddate = $date;
    }
  }
}

sub check_revline {

  my $line = shift;

  my ($date,$who,$comment);

  if( $line =~ /^[c\!]\s+revised\s+([\d\.]+)\s+by\s+(\S+)\s+(.*)$/i ) {
    $date = $1;
    $who = $2;
    $comment = $3;
  } elsif( $line =~ /^[c\!]\s+revised\s+on\s+([\d\.]+)\s+by\s+(\S+)\s+(.*)$/i ) {
    $date = $1;
    $who = $2;
    $comment = $3;
  } elsif( $line =~ /^[c\!]\s+revised\s+([\d\.]+)\s+(\S+)\s+(.*)$/i ) {
    $date = $1;
    $who = $2;
    $comment = $3;
  } elsif( $line =~ /^[c\!]\s+([\d\.]+)\s+(\S+)\s+(.*)$/i ) {
    $date = $1;
    $who = $2;
    $comment = $3;
  } else {
    #die "*** Cannot parse: $line\n";
    print STDERR "*** Cannot parse: $line\n";
  }

  return ($date,$who,$comment);
}






