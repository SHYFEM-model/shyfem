#!/usr/bin/perl -s

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

# use histo;
#
# my histo = new histo;
# $histo->init(1,2,3,4,5,6);	#set bins
# $histo->insert(4);
#
# my $nbins = $histo->get_nbins();
# my $bins  = $histo->get_bins();
# my $count = $histo->get_count();
#
# $histo->show("results");
# $histo->info("test");

#------------------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use strict;

package histo;

#-------------------------------------------------------------

sub new {

    my $self;

    $self =     {
			 nbins		=>	0
                         ,bins          =>      []
                         ,count         =>      []
		};

    bless $self;
    return $self;
}

#-------------------------------------------------------------

sub init {

  my ($self,@bins) = @_;

  my $n = scalar @bins;
  $self->{nbins} = $n;
  $self->{bins} = \@bins;
  my @count = ();
  $count[$n] = 0;
  @count = (0) x @count;
  $self->{count} = \@count;
}

sub insert {

  my ($self,$val) = @_;
  
  my $i = 0;
  for ( $i = 0; $i < $self->{nbins}; $i++ ) {
    last if( $self->{bins}->[$i] >= $val );
  }
  $self->{count}->[$i]++;
  #print STDERR "inserted $val at $i\n";
}

sub get_count {

  my ($self) = @_;

  return $self->{count};
}

sub get_nbins {

  my ($self) = @_;

  return $self->{nbins};
}

sub get_bins {

  my ($self) = @_;

  return $self->{bins};
}

sub show {

  my ($self,$text) = @_;

  $text = "" unless $text;

  print "show histogram: $text $self->{nbins}\n";

  my $oldlim = $self->{bins}->[0];
  print "x <= $oldlim: $self->{count}->[0]\n";
  for ( my $i = 1; $i < $self->{nbins}; $i++ ) {
    my $lim = $self->{bins}->[$i];
    print "$oldlim < x <= $lim: $self->{count}->[$i]\n";
    $oldlim = $lim;
  }
  my $lim = $self->{bins}->[-1];
  print "x > $lim: $self->{count}->[-1]\n";

}

sub info {

  my ($self,$text) = @_;
  
  $text = "" unless $text;

  print "info histogram: $text\n";
  print "nbin = $self->{nbins}\n";
  my $line = join(" ",@{$self->{bins}});
  print "bins = $line\n";
}

#-------------------------------------------------------------

sub histo_test {

  my $histo = new histo;
  $histo->init(1,2,3,4,5,6);
  $histo->info("test");
  $histo->insert(0);
  $histo->insert(1);
  $histo->insert(9);
  $histo->insert(5);
  $histo->insert(6);
  $histo->insert(4);
  $histo->show("results");
}

histo_test() if $0 =~ /histo.pm$/;

#-------------------------------------------------------------
1;
#-------------------------------------------------------------

