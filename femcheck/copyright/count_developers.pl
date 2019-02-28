#!/usr/bin/perl -ws

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use strict;

$::count = 0 unless $::count;
$::copy = 0 unless $::copy;
$::subst = 0 unless $::subst;
$::old = 0 unless $::old;
$::file = $ARGV[0];
%::devyear = ();
%::devcount = ();
%::names = ();

$::copystart = '!    Copyright (C)';

make_names();

count() if $::count;

while(<>) {

  chomp;
  if( $::old ) {
    check_old_revision();
  } else {
    check_revision();
  }
}

make_copyright() if $::copy;

#--------------------------------------------------------------

sub make_copyright
{
  my %hash = %::devcount;
  my @keys = sort { $hash{$b} <=> $hash{$a} } keys %hash;

  foreach my $key (@keys) {
     my $count = $::devcount{$key};
     my $years = consolidate_years($::devyear{$key});
     my $name = $::names{$key};
     #print "$key  $count   $years\n"; 
     print "$::copystart  $years  $name\n"; 
  }
}

sub consolidate_years {
  my $years = shift;

  $years =~ s/,$//;
  my @years = split(/,/,$years);
  my @new = ();

  my $oldyear = shift(@years);
  push(@new,$oldyear);
  foreach my $year (@years) {
    next if $year == $oldyear;
    push(@new,$year);
    $oldyear = $year;
  }

  my $line = join(",",@new);
  #print "$line\n";

  $line = "";
  my $start = shift(@new);
  my $end = $start;
  foreach my $year (@years) {
    if( $year - $end > 1 ) {
      if( $start == $end ) {
        $line .= "$end,";
      } else {
        $line .= "$start-$end,";
      }
      $start = $year;
      $end = $year;
    } else {
      $end = $year;
    }
  }
  if( $start == $end ) {
    $line .= "$end,";
  } else {
    $line .= "$start-$end,";
  }
  $line =~ s/,$//;

  return $line;
}

#--------------------------------------------------------------

sub check_old_revision {

  if( /^[cC!].*(\d{2}\.\d{2}\.\d{2})\s+/ ) {
    print "$_\n";
  }
}

sub check_revision {

  if( /^[cC!]\s+(\d{2}\.\d{2}\.\d{4})\s+(\S+)\s+/ ) {
    my $dev = $2;
    my $date = $1;
    my $year = $date;
    $year =~ s/^.*\.//;
    my @devs = split(/\&/,$dev);
    foreach my $d (@devs) {
      print "$d   $date  $::file\n" unless $::copy;
      $::devyear{$d} .= "$year,";
      $::devcount{$d}++;
    }
  }
}

sub count
{
  my %dev = ();

  while(<>) {
    my @f = split;
    my $d = $f[0];
    $dev{$d}++;
  }

  foreach my $d (sort keys %dev) {
    my $c = $dev{$d};
    my $name = $::names{$d};
    print "$d  $c  $name\n";
  }

  exit 0;
}

#--------------------------------------------------------------

sub make_names {

  %::names = (
         'ggu' => 'Georg Umgiesser'
        ,'aac' => 'Andrea Cucco'
        ,'aar' => 'Aaron Roland'
        ,'ccf' => 'Christian Ferrarin'
        ,'cpb' => ''
        ,'dbf' => 'Debora Bellafiore'
        ,'dmk' => 'Donata Melaku Canu'
        ,'erp' => 'Erik Pascolo'
        ,'fdp' => 'Francesca De Pascalis'
        ,'isa' => 'Isabella Scroccaro'
        ,'ivn' => 'Ivan Federico'
        ,'laa' => 'Leslie Aveytua'
        ,'lcz' => 'Lucia Zampato'
        ,'mbj' => 'Marco Bajo'
        ,'mcg' => 'Michol Ghezzo'
        ,'pzy' => 'Petras Zemlys'
        ,'wmk' => 'William McKiver'
    );
}

