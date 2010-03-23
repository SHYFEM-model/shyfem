#!/usr/bin/perl -s

# extracts entry from revision history on a specified date
#	and formats it to be included in gridlog
# can give -date on command line
# example:      "gridd -date 20-Oct-97 *.c"
# command option -write enables write of date to terminal

$line="";
$found=0;

if( $date ) {
  $compdate = shift;
} else {
  &usage;
}

while(<>) {
  if( /$compdate/ ) {
	&elabline($_);
	$found=1;
  } elsif( $found == 1 ) {
	if( /\d\d-[A-Z][a-z]{2}-\d\d/ ) {	#another date (matches any)
		$found=0;
	} elsif( /\s\*\s*\*/ ) {		#empty line
		$found=0;
	} else {				#continuation line
		&elabconti($_);
	}
  }
}

############ subroutines ##############

# print line

sub elabline {
  local($line,@dummy)=@_;
  if( $line =~ s/^\s*\*\s*(.*$compdate):// ) {
	$date = $1;
  }
  $line =~ s/^\s*//;
  $line =~ s/\s*\*\s*$/\n/;
  if( $write ) {
    print "- $ARGV: ($date)\n\t\t$line";
  } else {
    print "- $ARGV:\t$line";
  }
}

# print continuation line

sub elabconti {
  local($line,@dummy)=@_;
  $line =~ s/^\s*\*\s*(\S)/\1/;
  $line =~ s/\s*\*\s*$/\n/;
  print "\t\t  $line";
}

sub usage {
  print "Usage: gridd -date \"date-string\" file(s)\n";
  exit;
}
