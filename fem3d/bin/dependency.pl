#!/usr/bin/perl -s
#
# -f77   -ifort   -gfortran
#
#------------------------------------------

$file = $ARGV;

while(<>) {

  chomp;

  if( $f77 ) {
		f77();
	} elsif( $ifort ) {
		ifort();
  } else {
    print STDERR "unknown compiler or not specified\n";
    exit 1
  }

}

@subs = keys %subs;

foreach $prog (@subs) {

  print "  $prog\n";
}

exit 0;

#------------------------------------------------------

sub ifort {

  if( /^.* undefined reference to \`(\w+)\'$/ ) {
    $subs{$1}++;
  } elsif( /more undefined references/ ) {
    ;
  } elsif( /: In function \`/ ) {
    ;
  } elsif( /This statement function has not been used/ ) {
    ;
  } else {
    print STDERR "Cannot process: $_\n";
    exit 1
  }
}

sub f77 {

  if( /^.* undefined reference to \`(\w+)\'$/ ) {
    $subs{$1}++;
  } elsif( /more undefined references/ ) {
    ;
  } elsif( /collect2/ ) {
    ;
  } else {
    print STDERR "Cannot process: $_\n";
    exit 1
  }
}

