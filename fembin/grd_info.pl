#!/usr/bin/perl -w -s

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use strict;

$::dim = 1 if $::dim;

#---------------------------------------------------- options -----------

while( my $file = shift @ARGV ) {

  my $grid = new grd;
  $grid->{verbose} = 0;

  $grid->readgrd($file);				#FEM grid

  my $elems = $grid->get_elems();
  my $ne = scalar values %$elems;
  my $nodes = $grid->get_nodes();
  my $nn = scalar values %$nodes;
  my $lines = $grid->get_lines();
  my $nl = scalar values %$lines;

  my $ng = compute_grade($grid);

  if( $::dim ) {
    print STDERR "nkn = $nn    nel = $ne    ngr = $ng\n";
  } else {
    print STDERR "$file:   nodes $nn    elems $ne    lines $nl    grade $ng\n";
  }

}

sub compute_grade {

  my $grid = shift;

  my $grades = {};

  my $elems = $grid->get_elems();

  foreach my $elem (values %$elems) {
    my $vert = $elem->{vert};
    my $nold = $vert->[-1];
    foreach my $n (@$vert) {
      insert_grade($grades,$n,$nold);	# inserts connections into hash
      $nold = $n;
    }
  }

  my $ng = 0;
  foreach my $node (keys %$grades) {
    my $nlist = $grades->{$node};
    my $n = compute_grade_size($nlist);
    $ng = $n if $n > $ng;
  }

  return $ng;
}
      
sub insert_grade {

  my ($rhash,$n1,$n2) = @_;

  my $ra1 = $rhash->{$n1};
  unless( $ra1 ) {
    $ra1 = $rhash->{$n1} = [];
  }
  push(@$ra1,$n2);
 
  my $ra2 = $rhash->{$n2};
  unless( $ra2 ) {
    $ra2 = $rhash->{$n2} = [];
  }
  push(@$ra2,$n1);
}

sub compute_grade_size {	# eliminates double entries and returns size

  my $ra = shift;

  my @a = sort @$ra;
  my $na = @a;

  my $n = 1;
  for ( my $i=1; $i<$na; $i++ ) {
    $n++ if $a[$i-1] != $a[$i];
  }

  return $n;
}






