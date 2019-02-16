#!/usr/bin/perl -w

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use strict;

my %count = ();
my %files = ();

while(<>) {

  chomp;

  s/^\s+//;

  if( /^include\s*[\'\"](\w+)/ ) {
    my $inc = $1;
    #print "$inc\n";
    $count{$inc}++;
    $files{$inc} .= "$ARGV ";
  }
}

my @f = sort { $count{$a} <=> $count{$b} } keys %count;
my $tc = 0;
my $tcdim = 0;
my $incs = 0;

foreach my $key (@f) {
  my $c = $count{$key};
  my $file = $files{$key};
  my $cf = count_files($file);
  my $line = get_dims($key);
  #print "$c  $cf  $key ($file)";
  print "$c  $cf  $key";
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

sub count_files
{
  my $file = shift;

  my %files = ();
  my @f = split(/\s+/,$file);
  foreach (@f) {
    $files{$_}++;
  }
  my @files = keys %files;
  my $nf = @files;

  return $nf;
}

