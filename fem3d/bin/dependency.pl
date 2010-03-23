#!/usr/bin/perl

$file = $ARGV;

while(<>) {

  chomp;

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

@subs = keys %subs;

foreach $prog (@subs) {

  print "  $prog\n";
}

exit 0;
