#!/usr/bin/perl

my $skel = shift;
my $make = shift;

#------------------------- read in SKEL file

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

foreach my $key (keys %vals) {
  my $val = $vals{$key};
  print STDERR "$key  ->  $val\n";
}

#------------------------- substitute in Makefile

open(MAKE,"<$make") or die "Cannot open file: $make\n";

while(<MAKE>) {
  chomp;
  if( /^\s*(\w+)\s*=\s*/ ) {
    my $word = $1;
    foreach my $key (keys %vals) {
      if( $key eq $word ) {
        my $val = $vals{$key};
        print STDERR "substituting $key  ->  $val\n";
	$_ = "$key = $val";
      }
    }
  }
  print "$_\n";
}

#------------------------- finish
