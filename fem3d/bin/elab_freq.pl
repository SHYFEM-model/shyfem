#!/usr/bin/perl -w

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use strict;

my $name = "";
my @subs = ();

%::subs = ();
%::nsubs = ();


while(<>) {

  chomp;

  if( /^(\w+)/ ) {
    elab_file($name,\@subs);
    $name = $1;
    @subs = ();
  } elsif( /^\s+(\w+)/ ) {
    push(@subs,$1);
  } else {
    die "cannot parse: $_\n";
  }
}

elab_file($name,\@subs);

#----------------

my @keys = sort by_val keys %::nsubs;
@keys = reverse @keys;

foreach my $key (@keys) {
  my $n = $::nsubs{$key};
  print "$key  $n\n";
  if( $n > 0 and $n < 5 ) {
    my $subs = $::subs{$key};
    my $m = @$subs;
    foreach my $sub (@$subs) {
      print "  $sub";
    }
    print "\n";
  }
}

#-------------------------------------------------

sub elab_file
{
  my ($name,$subs) = @_;

  return unless $name;

  my @f = @$subs;
  my $n = @f;

  $::subs{$name} = \@f;
  $::nsubs{$name} = $n;
}

#-------------------------------------------------

sub by_val
{
  $::nsubs{$a} <=> $::nsubs{$b};
}
