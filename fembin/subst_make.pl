#!/usr/bin/perl -s
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# substitutes keywords and values in Makefile
#
# keywords and values can be given either in a file or on the command line
#
# first option:    subst_make.pl keyword-file makefile
#
# keyword-file:
#
# MACRO1 = VALUE1
# MACRO2 = VALUE2
# etc..
#
# second option:    
#
# subst_make.pl "MACRO1=VALUE1 MACRO2=VALUE2" makefile
#
# (no spaces around =)
#
#-------------------------------------------------------

$quiet = 1 if $quiet;	# do not be verbose
$first = 1 if $first;	# only change macro if it starts in first column

my $skel = shift;
my $make = shift;

my %vals = ();

#------------------------- read in SKEL file

if( -f $skel ) {		# read file
  %vals = read_skel($skel);
} else {			# macros on line
  %vals = parse_macros($skel);
}

foreach my $key (keys %vals) {
  my $val = $vals{$key};
  print STDERR "$key  ->  $val\n" unless $quiet;
}

#------------------------- substitute in Makefile

open(MAKE,"<$make") or die "*** Cannot open file: $make\n";

while(<MAKE>) {
  chomp;

  my $word = "";
  if( $first ) {		# must start in first column
    if( /^(\w+)\s*=\s*/ ) {
      $word = $1;
    }
  } else {
    if( /^\s*(\w+)\s*=\s*/ ) {
      $word = $1;
    }
  }

  if( $word ) {
    foreach my $key (keys %vals) {
      if( $key eq $word ) {
        my $val = $vals{$key};
        print STDERR "substituting $key  ->  $val\n" unless $quiet;
	$_ = "$key = $val";
      }
    }
  }
  print "$_\n";
}

#------------------------- subroutines

sub parse_macros {

  my $line = shift;

  my %vals = ();

  my @f = split(/\s+/,$line);

  foreach $macro (@f) {
    if( $macro =~ /^(\w+)=(\w+)$/ ) {
      my $key = $1;
      my $val = $2;
      $vals{$key} = $val;
    } else {
      die "*** Cannot parse macro: $macro\n";
    }
  }
  return %vals;
}

sub read_skel {

  my $skel = shift;

  my %vals = ();

  open(SKEL,"<$skel") or die "Cannot open file: $skel\n";

  while(<SKEL>) {
    chomp;
    if( /^\s*(\w+)\s*=\s*(.*)$/ ) {
      my $key = $1;
      my $val = $2;
      $vals{$key} = $val;
    } elsif( /^\s*$/ ) {		#skip empty line
    } else {
      die "Cannot parse line\n";
    }
  }

  close(SKEL);

  return %vals;
}

