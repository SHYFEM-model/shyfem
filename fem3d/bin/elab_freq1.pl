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

  if( /^-------/ ) {
    elab_file($name,\@subs);
    $name = <>; <>;
    $name = strip($name);
    @subs = ();
    #print "new file: $name\n"
  } elsif( /undefined reference to \`(\w+)_\'/ ) {
    push(@subs,$1);
    #print "  new sub: $1\n"
  }
}

elab_file($name,\@subs);

#----------------

exit;

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

  my %hash = ();

  foreach my $sub (@$subs) {
    $hash{$sub}++;
  }

  print "$name -------------------------\n";
  foreach my $key (sort keys %hash) {
    my $count = $hash{$key};
    print "$count    $key\n";
  }
}

#-------------------------------------------------

sub strip
{
  my $name = shift;

  chomp($name);
  $name =~ s/^ +//;
  $name =~ s/ +$//;

  return $name;
}

#-------------------------------------------------

sub by_val
{
  $::nsubs{$a} <=> $::nsubs{$b};
}
