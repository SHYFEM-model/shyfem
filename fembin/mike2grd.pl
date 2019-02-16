#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# converts from MIKE mesh format to GRD
#
#-------------------------------------------------------

#-------------------------------------------------------
# reading nodes
#-------------------------------------------------------

$_ = <>;
chomp;
my ($nkn,$proj) = split;

print STDERR "nkn,proj: $nkn $proj\n";

print"\n";
print"0 file transformed from MIKE mesh format\n";
print"\n";
print"0 projection: $proj\n";
print"\n";

for(my $i=1; $i<=$nkn; $i++) {
  $_ = <>;
  chomp;
  my ($k,$x,$y,$h,$t) = split;
  if( $i != $k ) {
    die "node number mismatch: $i $k\n";
  }
  my $d = -$h;
  print "1 $k $t $x $y $d\n";
}

#-------------------------------------------------------
# reading elements
#-------------------------------------------------------

$_ = <>;
chomp;
my ($nel,$nvert,$unknown) = split;

print STDERR "nel,nvert,unknown: $nel $nvert $unknown\n";
die "Expecting 3 vertices (triangles): $nvert\n" if $nvert != 3;
die "Expecting number 21: $unknown\n" if $unknown != 21;

for(my $i=1; $i<=$nel; $i++) {
  $_ = <>;
  chomp;
  my ($e,$v1,$v2,$v3) = split;
  if( $i != $e ) {
    die "element number mismatch: $i $e\n";
  }
  print "2 $e 0 3 $v1 $v2 $v3\n";
}

#-------------------------------------------------------
# final message
#-------------------------------------------------------

print STDERR "total number of nodes: $nkn\n";
print STDERR "total number of elements: $nel\n";
print STDERR "type of projection: $proj\n";
print STDERR "mesh file has been successfully converted...\n";

#-------------------------------------------------------
# end of routine
#-------------------------------------------------------

