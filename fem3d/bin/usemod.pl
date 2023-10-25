#!/usr/bin/perl -ws
#
# looks for files that need extra modules (not defined in the file)
#
#--------------------------------------------------------

use strict;

$::info = 0 unless $::info;
$::used = 0 unless $::used;
$::defined = 0 unless $::defined;
$::summary = 0 unless $::summary;

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

  if( scalar @list and not $::summary ) {
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

