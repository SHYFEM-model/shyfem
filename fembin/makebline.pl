#!/usr/bin/perl -w
#
# old version -> use bline.pl

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

  my $elems = $grid->get_elems();

  foreach my $elem (values %$elems) {
    my $verts = $elem->{vert};

    my $oldvert = $$verts[-1];
    foreach my $newvert (@$verts) {
      my $key = "$newvert:$oldvert";
      if( $keys{$key} ) {
        delete $keys{$key};
      } else {
        $key = "$oldvert:$newvert";
        $keys{$key} = 1;
      }
      $oldvert = $newvert
    }

    $grid->delete_elem($elem);			#here we delete the element
  }

  my $n = 0;
  foreach my $item (keys %keys) {
    my ($n1,$n2) = split(":",$item);
    $nodes{$n1} = $n2;
    #print STDERR "$item -> $n1  $n2\n";
    $n++;
  }
  print STDERR "Boundary nodes found: $n\n";

  my $nline = 0;
  foreach my $n1 (keys %nodes) {
    my $n2 = $nodes{$n1};
    if( $n2 ) {
      my $linenodes = make_line($n2,\%nodes);
      $nline++;
      my $nnodes = @$linenodes;
      #print STDERR "line $nline  with  $nnodes  nodes\n";

      my $line = $grid->make_line(0,0,0,@$linenodes);	#here we insert line
      $grid->close_line($line);
    }
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


