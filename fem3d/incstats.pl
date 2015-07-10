#!/usr/bin/perl -w

use strict;

my %count = ();

while(<>) {

  chomp;

  s/^\s+//;

  if( /^include\s*[\'\"](\w+)/ ) {
    my $inc = $1;
    #print "$inc\n";
    $count{$inc}++;
  }
}

my @f = sort { $count{$a} <=> $count{$b} } keys %count;
my $tc = 0;
my $tcdim = 0;
my $incs = 0;

foreach my $key (@f) {
  my $c = $count{$key};
  my $line = get_dims($key);
  print "$c   $key";
  if( $line ) {
    print "    ($line)";
    $tcdim++;
  }
  print "\n";
  $tc += $c;
  $incs++;
}

print "total: $incs   ($tcdim)   $tc\n";

#-------------------------------------------

sub get_dims {

  my $name = shift;

  my $file = $name . ".h";
  my %hash = ();

  unless ( open(FILE,"<$file") ) {
    print STDERR "*** cannot open file $file\n";
    return "";
  }

  #print STDERR "checking $file...\n";

  while( <FILE> ) {
    chomp;
    while( /^.*?(\w+dim\b)(.*)/ ) {
      $_ = $2;
      next if $1 eq "sedim";
      next if $1 eq "Sedim";
      next if $1 eq "bsedim";
      $hash{$1}++;
    }
  }

  close(FILE);

  my @keys = keys %hash;
  my $line = join(" ",@keys);

  return $line;
}

