#!/usr/bin/perl
#
# checks if file has dos line ending
#
#------------------------------------------------

my $file = "";
my $c = 0;

while(<>) {
  
  if( /\r\n$/ ) {
    if( $file ne $ARGV ) {
      print STDERR "$file $c\n" if $c;
      $c = 0;
      $file = $ARGV;
    }
    $c++;
  }
}

print STDERR "$file $c\n" if $c;

