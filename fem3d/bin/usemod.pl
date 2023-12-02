#!/usr/bin/perl -ws
#
# looks for files that need extra modules (not defined in the file)
#
#--------------------------------------------------------

use strict;

$::info = 0 unless $::info;
$::used = 0 unless $::used;
$::defined = 0 unless $::defined;
$::autosufficient = 0 unless $::autosufficient;
$::summary = 0 unless $::summary;
$::help = 0 unless $::help;
$::h = 0 unless $::h;
$::help = 1 if $::h;

if( $::help ) {
  print "Usage: usemod.pl [-help|-h] [-options] file\n";
  print "  options:\n";
  print "    -info           shows info on modules use and definition\n";
  print "    -used           shows used modules\n";
  print "    -defined        shows defined modules\n";
  print "    -autosufficient shows modules that are autosufficient\n";
  print "    -summary        shows summary of used and defined modules\n";
  print "    -h|-help        this help screen\n";
  exit 0;
} elsif( !@ARGV ) {
  print "Usage: usemod.pl [-help|-h] [-options] file\n";
  exit 0;
}

my $file = "";
%::defined = ();
%::used = ();
@::defined = ();	# for summary
@::used = ();

while(<>) {

  if( $file ne $ARGV ) {
    print_info($file);
    $file = $ARGV;
  }

  chomp;

  if( /^\s*module\s+(\w+)/ ) {		# definition of module
    $::defined{$1}++;
  } elsif( /^\s*use\s+(\w+)/ ) {		# usage of module
    $::used{$1}++;
  }

}

print_info($file);
print_summary();

#--------------------------------------------------------

sub print_info
{
  my $file = shift;

  return unless $file;

  if( scalar %::used ) {
    if( $::info or $::used ) {
      print "modules used in file $file:\n";
      foreach my $m (sort keys %::used) {
        print "  $m\n";
      }
    }
    push(@::used,keys %::used);
  }

  if( scalar %::defined ) {
    if( $::info or $::defined ) {
      print "modules defined in file $file:\n";
      foreach my $m (sort keys %::defined) {
        print "  $m\n";
      }
    }
    push(@::defined,keys %::defined);
  }

  my @list = ();
  foreach my $m (keys %::used) {
    push(@list,$m) unless $::defined{$m};
  }

  if( not scalar @list and $::autosufficient ) {
    print "autosufficient usage of modules in file $file\n";
  }

  if( scalar @list and not $::summary and not $::autosufficient ) {
    print "modules used but not defined in file $file:\n";
    foreach my $m (sort @list) {
      print "  $m\n";
    }
  }

  %::defined = ();
  %::used = ();
}

#--------------------------------------------------------

sub print_summary
{
  return unless $::summary;

  if( scalar @::used or scalar @::defined ) {
    #print "======== summary =========\n";
  }

  if( scalar @::used ) {
    my %list = ();
    foreach my $m (@::used) {
      $list{$m}++;
    }
    print "summary of modules used in all files:\n";
    foreach my $m (sort keys %list) {
      print "  $m\n";
    }
  }

  if( scalar @::defined ) {
    print "summary of modules defined in all files:\n";
    foreach my $m (sort @::defined) {
      print "  $m\n";
    }
  }
}

#--------------------------------------------------------

