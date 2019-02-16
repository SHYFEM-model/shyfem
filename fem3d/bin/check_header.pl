#!/usr/bin/perl -w

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use strict;

%::first = ();
%::keywords = ();

foreach my $file (@ARGV) {
  treat_files($file);
}

write_hash("first statement",\%::first);
write_hash("key words",\%::keywords);

#---------------------------------------------------------

sub treat_files {

  my $file = shift;

  open(FILE,"<$file");
  my $text = get_header($file);
  my $rhash = get_keywords($text);
  my $rlist = get_revision($text);
  print "no keywords in: $file\n" unless $rhash;
  insert_hash($rhash);
  close(FILE);
}

#---------------------------------------------------------

sub get_revision {

  my $text = shift;

  my $n = 0;
  my %hash = ();

  foreach my $line (@$text) {
    if( $line =~ /^[!cC]\s+(.*)\s+:\s*$/ ) {
sub get_keywords {

  my $text = shift;

  my $n = 0;
  my %hash = ();

  foreach my $line (@$text) {
    #print "$line\n";
    if( $line =~ /^[!cC]\s+(.*)\s+:\s*$/ ) {
      my $key = $1;
      $hash{$key}++;
      $n++;
    }
  }

  return \%hash;
}

sub get_header {

  my $file = shift;

  my @header = ();

  while(<FILE>) {

    chomp;

    if( /^\s+(\w+)/ ) {
      my $name = lc($1);
      $::first{$name}++;
      #print STDERR "******* $file\n" if $name eq "include";
      #print STDERR "******* $file\n" if $name eq "integer";
      last;
    }

    push(@header,$_);
  }

  return \@header;
}

#---------------------------------------------------------

sub insert_hash {

  my $rhash = shift;

  foreach my $key (keys %$rhash) {
    $::keywords{$key}++;
  }
}

sub write_hash {

  my ($text,$rhash) = @_;

  print STDERR "$text:\n";
  foreach my $key (keys %$rhash) {
    my $val = $rhash->{$key};
    print STDERR "   $val   $key\n";
  }
}

#---------------------------------------------------------

