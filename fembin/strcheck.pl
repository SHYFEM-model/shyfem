#!/usr/bin/perl -w

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use date;
use strict;
 
my $date = new date;

my $debug = 1;

my $err = 0;
my $warn = 0;
my ($file,$tanf,$tend,$nrecs);

my $itanf = shift;
my $itend = shift;
my $date0 = shift;

die "itanf undefined\n" unless defined $itanf;
die "itend undefined\n" unless defined $itend;
die "date0 undefined\n" unless defined $date0;

$date->init_date($date0);
$itanf = to_it($date,$itanf);
$itend = to_it($date,$itend);

print "STR params: $itanf $itend $date0\n" if $debug;

while(<>) {
  chomp;

  print "debug: $_\n" if $debug;

  if( /file name:\s*(\S+)/ ) { $file = $1; }
  if( /nrecs:\s*(\S+)/ ) { $nrecs = $1; }
  if( /start time:\s*(\S+)/ ) { $tanf = $1; }
  if( /end time:\s*(\S+)/ ) { $tend = $1; }
  if( /^\s*\*\*\* / ) {
    $err++;
    print STDERR "$_\n";
  } elsif( /^\s*\* / ) {
    $warn++;
    print STDERR "$_\n";
  }
}

die "tanf undefined\n" unless defined $tanf;
die "tend undefined\n" unless defined $tend;
die "nrecs undefined\n" unless defined $nrecs;
die "file undefined\n" unless defined $file;

print "file params: $tanf $tend $nrecs $file\n" if $debug;

if( $nrecs > 1 ) {
  $err += check_time_params($file,$date,$itanf,$itend,$tanf,$tend);
}

if( $err ) {
  print STDERR "*** there were errors in file $file\n";
} elsif( $warn ) {
  print STDERR "* there were warnings in file $file\n";
}

if( $err ) {
  exit 1;
} elsif( $warn ) {
  exit 2;
} else {
  exit 0
}

#---------------------------------------------------------

sub check_time_params {

  my ($file,$date,$itanf,$itend,$tanf,$tend) = @_;

  my $err = 0;
  my $line = "";

  $itanf = to_abs($date,$itanf);
  $itend = to_abs($date,$itend);
  $tanf = to_abs($date,$tanf);
  $tend = to_abs($date,$tend);

  if( $itanf < $tanf ) {
    print STDERR "*** start of simulation before start of file\n";
    $err = 1;
  } elsif( $itend > $tend ) {
    print STDERR "*** end of simulation after end of file\n";
    $err = 1;
  }
  if( $err ) {
    print STDERR "file: $file\n";
    $line = $date->format_abs($itanf);
    print STDERR "start of simulation: $itanf    $line\n";
    $line = $date->format_abs($tanf);
    print STDERR "start of file      : $tanf    $line\n";
    $line = $date->format_abs($itend);
    print STDERR "end of simulation  : $itend    $line\n";
    $line = $date->format_abs($tend);
    print STDERR "end of file        : $tend    $line\n";
  }
}

#---------------------------------------------------------

sub to_it {

  my ($date,$time) = @_;

  if( $time =~ /^\'(.*)\'$/ ) {		# is a string
    my $string = $1;
    my @f = $date->unformat_time_date($string);
    $time = $date->convert_to_it(@f);
  }

  return $time;
}

sub to_abs {

  my ($date,$rel) = @_;

  my @f = $date->convert_from_it($rel);
  my $abs = $date->convert_to_abs(@f);

  return $abs;
}

#---------------------------------------------------------

