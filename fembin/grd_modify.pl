#!/usr/bin/perl -w -s
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# version = 2.4
#
# 07.10.2010	ggu	act on all items if no line is given
# 07.10.2010	ggu	translate nodes
# 10.02.2012	ggu	-depth_invert
# 16.02.2018	ggu	-unset_depth
# 23.04.2019	ggu	-compress
# 25.06.2019	ggu	-nodes
# 27.04.2021	ggu	-depth_diff
# 29.04.2021	ggu	-join_lines and -min_area
# 30.04.2021	ggu	-clean and -unused
# 15.06.2022	ggu	new -write for elems and lines, set must_invert
#
#----------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use grdline;
use strict;

#---------------------------------------------------- options -----------

$::help = 1 if $::help;
$::h = 1 if $::h;
$::invert = 0 unless $::invert;
$::exclude = 0 unless $::exclude;
$::preserve = 0 unless $::preserve;
$::compress = 0 unless $::compress;
$::depth_invert = 0 unless $::depth_invert;
$::write = 0 unless $::write;
$::verbose = 0 unless $::verbose;
$::nodes = "" unless $::nodes;
$::join_lines = "" unless $::join_lines;
$::min_area = 0 unless $::min_area;
$::min_dist = 0 unless $::min_dist;
$::clean = "" unless $::clean;
$::unused = "" unless $::unused;

if( $::n or $::e or $::l ) {	#explicitly given -> only on these items
  $::n = 0 unless $::n;
  $::e = 0 unless $::e;
  $::l = 0 unless $::l;
} else {			#default -> on all items
  $::n = 1;
  $::e = 1;
  $::l = 1;
}

$::flag = -99999;
$::depth_flag = -999;

$::print = 0 unless $::print;
$::delete = 0 unless $::delete;
$::type = $::flag unless defined($::type);
$::depth = $::flag unless defined($::depth);
$::depth_diff = $::flag unless defined($::depth_diff);
$::trans = "" unless defined($::trans);
$::latlon = "" unless defined($::latlon);

$::unused = 1 if $::min_area;

if( $::nodes ) {
  @::nodes = split(",",$::nodes);
}

if( $::trans ) {
  ($::trans_x,$::trans_y) = split(",",$::trans);
  print STDERR "Translating with $::trans_x $::trans_y\n";
}

$::inside = 1;		#if true delete inside
$::outside = 0;		#if true delete outside

if( $::invert ) {
  $::outside = 1;
  $::inside = 0;
}

if( $::h or $::help ) {
  FullUsage();
  exit 0;
} elsif( not $ARGV[0] ) {
  Usage();
  exit 1;
} elsif( not $ARGV[1] ) {
  #Usage();
  #exit 1;
  print STDERR "No line given -> applying selection to all items\n";
}

print STDERR "in: $::inside  out: $::outside  exclude: $::exclude\n";
print STDERR "print: $::print  delete: $::delete\n";
print STDERR "depth: $::depth  type: $::type\n";
print STDERR "depth_diff: $::depth_diff\n";
print STDERR "n: $::n  e: $::e  l: $::l\n";

#------------------------------------------------- read files -----------

my $grid = new grd;
my $file = $ARGV[0];
$grid->readgrd($file);				#FEM grid

my $gline = new grd;
my $gfile = $ARGV[1];
if( $gfile ) {
  $gline->readgrd($gfile);			#grid with lines
} else {
  make_dummy_line($grid,$gline);		# dummy line around all nodes
}

#-------------------------------------------------- main ----------------

$grid->set_latlon($::latlon);
make_end_lines($grid) if $::join_lines;

my $lines = $gline->get_lines();

foreach my $line (values %$lines) {
  my $grline = make_grdline($gline,$line);	#sets up new datastructure
  flag_nodes($grid,$grline);			#sets up $::flags
  #print_flaged_nodes();

  loop_on_elements($grid,\&modify_element) if $::e;
  loop_on_lines($grid,\&modify_line) if $::l;
  loop_on_nodes($grid,\&modify_node) if $::n;
}

$grid->delete_degenerate() if $::clean;
$grid->delete_unused() if $::unused;

make_compress($grid) if $::compress;
$grid->set_preserve_order($::preserve);
$grid->writegrd("modify.grd");

#------------------------------------------------------------------------

sub FullUsage {
  print STDERR "                                    \n";
  Usage();
  print STDERR "                                    \n";
  print STDERR "  Selects and modifies items of grid contained in line\n";
  print STDERR "                                    \n";
  print STDERR "  If no line is given acts on all items of grid\n";
  print STDERR "                                    \n";
  print STDERR "  -h|-help      this help screen\n";
  print STDERR "                                    \n";
  print STDERR "  {-n|-e|-l}    modify nodes,elements,lines (default: all)\n";
  print STDERR "  -verbose      write more messages to terminal\n";
  print STDERR "  -print        print selected items\n";
  print STDERR "  -write        write element or line numbers to terminal\n";
  print STDERR "  -type=type    set type of selected items to type\n";
  print STDERR "  -depth=depth  set depth of selected items to depth\n";
  print STDERR "  -depth_diff=d set depth of selected items to depth+d\n";
  print STDERR "  -trans=dx,dy  translate selected nodes by dx/dy\n";
  print STDERR "  -nodes=list   extract info about given nodes in list\n";
  print STDERR "  -delete       delete selected items\n";
  print STDERR "                                    \n";
  print STDERR "  -invert       invert selection - modify outer area\n";
  print STDERR "  -exclude      exclude border line items from selection\n";
  print STDERR "  -preserve     preserves order of items\n";
  print STDERR "  -depth_invert inverts depth values (neg to pos etc.)\n";
  print STDERR "  -unset_depth  deletes depth values\n";
  print STDERR "  -compress     compresses node and element numbers\n";
  print STDERR "  -min_area=a   deletes islands with area < a\n";
  print STDERR "  -join_lines   join lines that end at same node\n";
  print STDERR "  -min_dist=d   joins lines with final node distance < d\n";
  print STDERR "  -unify_node=d unifies nodes with distance < d\n";
  print STDERR "  -latlon       treat coordinates as lat/lon\n";
  print STDERR "  -clean        cleans grid from degenerate elements/lines\n";
  print STDERR "  -unused       deletes unused nodes\n";
  print STDERR "                                    \n";
  print STDERR "islands are closed lines\n";
}

sub Usage {
  print STDERR "Usage: grd_modify.pl [-h|-help] [-options] grid [line]\n";
}

#----------------------------------------------------------

sub must_handle_item
{
  my ($item) = @_;

  return 0 unless $item;

  my $vert = {};

  if( exists($item->{vert}) ) {		#is elem or line
    $vert = $item->{vert};
  } else {				#is node - insert only one value
    my @vert = ($item->{number}); $vert = \@vert;
  }

  my $nvert = @$vert;

  my $n = 0;
  foreach my $node (@$vert) {
    $n++ if $::flags{$node};
  }

  if( $::inside ) {
      if( $::exclude ) {
        return 1 if( $n == $nvert );
      } else {
        return 1 if( $n > 0 );
      }
  } else {
      if( $::exclude ) {
        return 1 if( $n == 0 );
      } else {
        return 1 if( $n < $nvert );
      }
  }

  return 0;
}

sub internal_check {	# checks if both approaches give same answer

  my ($item) = @_;

  my $number = $item->{number};
  my $f = $::flags{$number};

  my $y1 = must_handle_item($item);
  my $y2 = ( $f and $::inside or not $f and not $::inside );
  $y2 = 1 if $y2;

  die "internal check error: $y1 $y2\n" if $y1 != $y2;
}

#-----------------------------------------------------------------

sub loop_on_nodes {

  my ($grid,$proc) = @_;

  $grid->make_used();
  my $nodes = $grid->get_nodes();

  foreach my $node (values %$nodes) {
    internal_check($node);
    if( must_handle_item($node) ) {
      &$proc($grid,$node);
    }
  }

  foreach my $nnumber ( @::nodes ) {
    my $node = $nodes->{$nnumber};
    my $x = $node->{x};
    my $y = $node->{y};
    my $h = $node->{h};
    $h = "" if not $grid->has_depth($node);
    print "1 $nnumber 3 $x $y $h\n";
  }
}

sub loop_on_elements {

  my ($grid,$proc) = @_;

  my $elems = $grid->get_elems();

  loop_on_items($grid,$proc,$elems);
}
 
sub loop_on_lines {

  my ($grid,$proc) = @_;

  my $lines = $grid->get_lines();

  loop_on_items($grid,$proc,$lines);
}

sub loop_on_items {

  my ($grid,$proc,$items) = @_;

  foreach my $item (values %$items) {
    if( must_handle_item($item) ) {
      &$proc($grid,$item);
    }
  }
}

#-----------------------------------------------------------------

sub modify_node {

  my ($grid,$node) = @_;

  my $number = $node->{number};

  if( $::print ) {
    print "modifying node $number\n";
  } elsif( $::type >= 0 ) {
    $node->{type} = $::type;
  } elsif( $::delete ) {
    $grid->delete_node($node) unless $node->{used};
  } elsif( $::depth != $::flag ) {
    $node->{h} = $::depth;
  } elsif( $::depth_diff != $::flag ) {
    $node->{h} += $::depth_diff if $grid->has_depth($node);
  } elsif( $::depth_invert ) {
    $node->{h} = -$node->{h};
  } elsif( $::unset_depth ) {
    $node->{h} = $::depth_flag;
  } elsif( $::trans ) {
    $node->{x} += $::trans_x;
    $node->{y} += $::trans_y;
  }
}

sub modify_element {

  my ($grid,$elem) = @_;

  my $number = $elem->{number};

  if( $::print ) {
    print "modifying element $number\n";
  } elsif( $::write ) {
    write_item($grid,$elem);
  } elsif( $::type != $::flag ) {
    $elem->{type} = $::type;
  } elsif( $::delete ) {
    $grid->delete_elem($elem);
  } elsif( $::depth != $::flag ) {
    $elem->{h} = $::depth;
  } elsif( $::depth_diff != $::flag ) {
    $elem->{h} += $::depth_diff if $grid->has_depth($elem);
  } elsif( $::depth_invert ) {
    $elem->{h} = -$elem->{h};
  } elsif( $::unset_depth ) {
    $elem->{h} = $::depth_flag;
  }
}

sub modify_line {

  my ($grid,$line) = @_;

  my $number = $line->{number};

  if( $::print ) {
    print "modifying line $number\n";
  } elsif( $::write ) {
    write_item($grid,$line);
  } elsif( $::type != $::flag ) {
    $line->{type} = $::type;
  } elsif( $::delete ) {
    $grid->delete_line($line);
  } elsif( $::depth != $::flag ) {
    $line->{h} = $::depth;
  } elsif( $::depth_diff != $::flag ) {
    $line->{h} += $::depth_diff if $grid->has_depth($line);
  } elsif( $::depth_invert ) {
    $line->{h} = -$line->{h};
  } elsif( $::unset_depth ) {
    $line->{h} = $::depth_flag;
  } elsif( $::join_lines ) {
    join_line($grid,$line);
  } elsif( $::min_area ) {
    delete_island($grid,$line);
  }
}

#-----------------------------------------------------------------

sub write_item {

  my ($grid,$item) = @_;

  return unless must_handle_item($item);

  my $vert = $item->{vert};

  foreach my $node (@$vert) {
    print "$node\n";
  }
}

#-----------------------------------------------------------------

sub delete_island {

  my ($grid,$litem) = @_;

  return unless must_handle_item($litem);
  return unless $grid->is_closed($litem);

  my $area = $grid->area($litem);
  if( $area < $::min_area ) {
    my $nl = $litem->{number};
    print STDERR "deleting island $nl with area $area\n";
    $grid->delete_line($litem);
  }
}

sub join_line {

  my ($grid,$litem) = @_;

  return unless must_handle_item($litem);
  return if $grid->is_closed($litem);

  my $node;
  my $litem2 = "";
  my $vert = $litem->{vert};
  if( $litem2 = find_other_line($litem,$vert->[0]) ) {
    $node = $vert->[0];
  } elsif( $litem2 = find_other_line($litem,$vert->[-1]) ) {
    $node = $vert->[-1];
  } elsif( $litem2 = find_close_line($grid,$litem,$vert->[0]) ) {
    $node = $vert->[0];
  } elsif( $litem2 = find_close_line($grid,$litem,$vert->[-1]) ) {
    $node = $vert->[-1];
  }
  return unless must_handle_item($litem2);
  adjust_end_lines($node,$litem,$litem2),
  #delete_end_lines($node);

  my $nl1 = $litem->{number};
  my $nl2 = $litem2->{number};

  print STDERR "joining lines $nl1 and $nl2\n";

  my $vert1 = $litem->{vert};
  my $vert2 = $litem2->{vert};
  my $nv1 = @$vert1;
  my $nv2 = @$vert2;
  $grid->connect_lines($litem2,$litem);	#first litem2 to avoid freed memory bug
  my $vert3 = $litem2->{vert};
  my $nv3 = @$vert3;
  my $ndiff = $nv3 - $nv1 - $nv2 + 1;
}
 
sub find_other_line {

  my ($litem,$node) = @_;

  return unless defined $::end_lines{$node};
  my $ends = $::end_lines{$node};
  my $n = @$ends;

  if( $n == 2 ) {
    if( $ends->[0] == $litem ) {
      return $ends->[1];
    } elsif( $ends->[1] eq $litem ) {
      return $ends->[0];
    } else {
      die "internal error joining...\n";
    }
  } elsif( $n > 3 ) {
    print STDERR "$n lines ending on node $node...\n";
  }
  return "";
}

sub find_close_line {

  my ($grid,$litem,$node) = @_;

  return "" if $::min_dist <= 0;	#only if positive distance

  return unless defined $::end_lines{$node};
  my $ends = $::end_lines{$node};
  my $n = @$ends;
  return unless $n == 1;

  foreach my $nend (keys %::end_lines ) {
    next if $node == $nend;
    $ends = $::end_lines{$nend};
    $n = @$ends;
    #print STDERR ".... $node  $::min_dist  $nend $n\n";
    next unless $n == 1;
    my $item = $ends->[0];
    next if $item == $litem;
    if( $grid->dist($node,$nend) < $::min_dist ) {
      my $vert = $item->{vert};
      if( $vert->[0] == $nend ) {
        $vert->[0] = $node;
      } elsif( $vert->[-1] == $nend ) {
        $vert->[-1] = $node;
      } else {
        print STDERR ".... $node $nend $::min_dist @$vert\n";
        die "internal error find_close_line\n";
      }
      delete_end_lines($nend);
      return $item;
    }
  }

  return "";
}

#-----------------------------------------------------------------

sub make_compress {

  my ($grid) = @_;

  print STDERR "compressing item numbers...\n";

  my $n;
  my $items;
  my %node_numbers = ();

  my $nodes = $grid->get_nodes();
  $n = 0;
  my %nnew = ();
  foreach my $node (values %$nodes) {
    my $number = $node->{number};
    $node->{number} = ++$n;
    $node_numbers{$number} = $n;
    $nnew{$n} = $node;
  }
  $grid->{nodes} = \%nnew;

  $items = $grid->get_elems();
  $n = 0;
  my %enew = ();
  foreach my $item (values %$items) {
    $item->{number} = ++$n;
    my $vert = $item->{vert};
    foreach my $node (@$vert) {
      $node = $node_numbers{$node};
    }
    $enew{$n} = $item;
  }
  $grid->{elems} = \%enew;

  $items = $grid->get_lines();
  $n = 0;
  my %lnew = ();
  foreach my $item (values %$items) {
    $item->{number} = ++$n;
    my $vert = $item->{vert};
    foreach my $node (@$vert) {
      $node = $node_numbers{$node};
    }
    $lnew{$n} = $item;
  }
  $grid->{lines} = \%lnew;
}

#-----------------------------------------------------------------

sub print_flaged_nodes {

  foreach my $number (keys %::flags) {
    print "$number\n";
  }
}

sub flag_nodes {

  my ($grid,$grline) = @_;

  my $nin = 0;
  my $nout = 0;

  %::flags = ();

  my $nodes = $grid->get_nodes();
  foreach my $node (values %$nodes) {
    my $number = $node->{number};
    my $x = $node->{x};
    my $y = $node->{y};
    if( $grline->in_line($x,$y) ) {
      $::flags{$number} = 1;
      $nin++;
    } else {
      $nout++;
    }
  }
  print STDERR "flag_nodes:  in: $nin  out: $nout\n" if $::verbose;
}

#-----------------------------------------------------------------

sub make_grdline {

  my ($gline,$line) = @_;

  my ($x,$y) = make_xy_array($gline,$line);
  my $grline = new grdline;
  $grline->{must_invert} = 1;
  $grline->set_line($x,$y);

  return $grline;
}

sub make_xy_array {

  my ($grid,$line) = @_;

  my @x = ();
  my @y = ();

  my $vert = $line->{vert};
  foreach my $nnode (@$vert) {
    my $node = $grid->get_node($nnode);
    my $x = $node->{x};
    my $y = $node->{y};
    push(@x,$x);
    push(@y,$y);
  }

  return (\@x,\@y);
}

#-----------------------------------------------------------------

sub make_dummy_line {

  my ($grid,$gline) = @_;

  my ($xmin,$ymin,$xmax,$ymax) = $grid->get_xy_minmax();

  print STDERR "xy_minmax: $xmin $ymin $xmax $ymax \n";

  $xmin -= 0.1*($xmax-$xmin);
  $xmax += 0.1*($xmax-$xmin);
  $ymin -= 0.1*($ymax-$ymin);
  $ymax += 0.1*($ymax-$ymin);

  $gline->make_node(1,0,$xmin,$ymin);
  $gline->make_node(2,0,$xmax,$ymin);
  $gline->make_node(3,0,$xmax,$ymax);
  $gline->make_node(4,0,$xmin,$ymax);

  $gline->make_line(1,0,0,1,2,3,4,1);
}

#-----------------------------------------------------------------

sub list2hash {

  my @a = @_;

  my %h = ();

  foreach (@a) {
    $h{$_} = 1;
  }

  return %h
}

#-----------------------------------------------------------------

sub make_end_lines {

  my $grid = shift;

  %::end_lines = ();

  my $litems = $grid->get_lines();

  foreach my $litem (values %$litems) {
    my $vert = $litem->{vert};
    my $nv = @$vert;
    next if $grid->is_closed($litem);
    insert_end_lines($vert->[0],$litem);
    insert_end_lines($vert->[-1],$litem);
  }
}

sub insert_end_lines {

  my ($n,$litem) = @_;

  if( not exists($::end_lines{$n}) ) {
    my @a = ();
    $::end_lines{$n} = \@a;
  }

  my $a = $::end_lines{$n};
  push(@$a,$litem);
}

sub delete_end_lines {

  my ($n) = @_;

  delete $::end_lines{$n};
}

sub adjust_end_lines {

  my ($node,$item1,$item2) = @_;		#item1 will be deleted

  delete $::end_lines{$node};

  my $vert = $item1->{vert};

  my $node2;
  if( $vert->[0] == $node ) {
    $node2 = $vert->[-1];
  } elsif( $vert->[-1] == $node ) {
    $node2 = $vert->[0];
  } else {
    die "internal error adjust_end_lines\n";
  }

  my $ends = $::end_lines{$node2};
  foreach my $end (@$ends) {
    $end = $item2 if $end == $item1; 
  }
}


#-----------------------------------------------------------------

#-----------------------------------------------------------------

