#!/usr/bin/perl -w

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use strict;

my $grid = new grd;
my $file = $ARGV[0];

$grid->readgrd($file);

$grid = make_bound_line($grid);

$grid->writegrd;

#######################################################################

sub make_bound_line
{
  my $grid = shift;
  my %keys = ();
  my %nodes = ();
  my %bnodes = ();

  my $rkeys = \%keys;

  my $elems = $grid->get_elems();

  foreach my $elem (values %$elems) {
    my $verts = $elem->{vert};

    my $oldvert = $$verts[-1];
    foreach my $newvert (@$verts) {
      if( not delete_keys($rkeys,$oldvert,$newvert) ) {
        insert_keys($rkeys,$oldvert,$newvert);
      }
      $oldvert = $newvert;
    }

    $grid->delete_elem($elem);			#here we delete the element
  }

  my $n = 0;
  my $error = 0;
  foreach my $n1 (keys %keys) {
    my $nodes = $keys{$n1};
    $nodes =~ s/^\s+//;
    $nodes =~ s/\s+$//;
    $nodes =~ s/\s+/ /;
    my @f = split(/\s+/,$nodes);
    my $nn = @f;
    if( $nn ne 2 and $nn ne 0 ) {
      print STDERR "more than two lines from node $n1 ($nn): $nodes\n";
      $error++;
    } elsif( $nn == 2 ) {
      my ($n2,$n3) = @f;
      if( $n2*$n3 >= 0 ) {
        print STDERR "error in structure: $n2 $n3\n";
      } else {
	if( $n2 > 0 ) {
          $bnodes{$n1} = $n2;
	} else {
          $bnodes{$n1} = $n3;
	}
      }
    }
    $n++;
  }
  if( $error ) {
    die "$error errors in boundary line found. Aborting.\n";
  }
  print STDERR "Boundary nodes found: $n\n";

  my $nline = 0;
  foreach my $n1 (keys %bnodes) {
    my $n2 = $bnodes{$n1};
    if( $n2 ) {
      my $linenodes = make_line($n2,\%bnodes);
      $nline++;
      my $nnodes = @$linenodes;
      #print STDERR "line $nline  with  $nnodes  nodes\n";

      my $line = $grid->make_line(0,0,0,@$linenodes);   #here we insert line
      $grid->close_line($line);
    }

    my $nodes = $keys{$n1};
  }
  print STDERR "Boundary lines found: $nline\n";

  $grid->delete_unused();			#delete unused nodes

  return $grid;
}

sub make_line
{
  my ($start,$nodes) = @_;

  my @line = ();
  my $old = $start;
  my $node;

  while( $node = $$nodes{$old} ) {
    push(@line,$node);
    $$nodes{$old} = 0;
    $old = $node;
  }

  if( $old != $start ) {
    die "error in line: $old  $start\n";
  }

  return \@line;
}

#----------------------------------------------------

sub insert_keys {

  my ($rkeys,$oldvert,$newvert) = @_;

  $$rkeys{$oldvert} .= " $newvert ";
  $$rkeys{$newvert} .= " -$oldvert ";
}

sub delete_keys {

  my ($rkeys,$oldvert,$newvert) = @_;

  my $found = 0;
  my $error = 0;

  if( defined $$rkeys{$oldvert} and $$rkeys{$oldvert} =~ / -$newvert / ) {
    $found = 1;
    $$rkeys{$oldvert} =~ s/ -$newvert //;
  }

  if( defined $$rkeys{$newvert} and $$rkeys{$newvert} =~ / $oldvert / ) {
    if( not $found ) {
      $error = 1;
    } else {
      $$rkeys{$newvert} =~ s/ $oldvert //;
    }
  } elsif( $found ) {
    $error = 1;
  }

  if( $error ) {
    die "error in keys structure\n";
  }

  return $found;
}

#----------------------------------------------------

