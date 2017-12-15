#!/usr/bin/perl
#
# compacts data structure routines into one file
#
#---------------------------------------------------------

while( my $file = shift(@ARGV) ) {

  print STDERR "$file\n";

  compress_file($file);
}

sub compress_file {

  my $file = shift;

  open(FILE,"<$file");

  print "!*********************************************************\n";
  print "!*********************************************************\n";
  print "!*********************************************************\n";
  print "! $file\n";
  print "!*********************************************************\n";
  print "!*********************************************************\n";
  print "!*********************************************************\n";

  while(<FILE>) {
    print;
    last if /^\s+end\s+module/;
  }
  $_ = <FILE>; print; print "\n";

  close(FILE);
}

