#!/usr/bin/perl
#
# splits path in single parts and compares with given value
# path is a : delimited string with no spaces
#
# Usage: split_path.pl compare path
#
#-------------------------------------------------

$debug = 0;

$compare = $ARGV[0];
$path = $ARGV[1];

print STDERR "looking for $compare in path $path\n" if $debug;

Exit(0) unless $compare;
Exit(0) unless $path;

@path = split(/:/,$path);

foreach $item (@path) {
  if( $item eq $compare ) {
    print STDERR "found $compare in path\n" if $debug;
    Exit(1);
  }
}

print STDERR "no $compare found in path\n" if $debug;
Exit(0);

#-------------------------------------------------

sub Exit()
{
  my $status = shift;

  print "$status\n";
  exit $status;
}

#-------------------------------------------------

