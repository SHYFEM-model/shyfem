#!/usr/bin/perl -w -s
#
# version = 2.1
#
# 07.10.2010	ggu	act on all items if no line is given
# 07.10.2010	ggu	translate nodes
# 10.02.2012	ggu	-depth_invert
# 16.02.2018	ggu	-unset_depth
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
$::depth_invert = 0 unless $::depth_invert;

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
$::trans = "" unless defined($::trans);

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

my $lines = $gline->get_lines();
my $flag = {};

foreach my $line (values %$lines) {
  my $grline = make_grdline($gline,$line);	#sets up new datastructure
  flag_nodes($grid,$grline,$flag);
  #print_flaged_nodes($flag);

  loop_on_elements($grid,$flag,\&modify_element) if $::e;
  loop_on_lines($grid,$flag,\&modify_line) if $::l;
  loop_on_nodes($grid,$flag,\&modify_node) if $::n;
}

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
  print STDERR "  -print        print selected items\n";
  print STDERR "  -type=type    set type of selected items to type\n";
  print STDERR "  -depth=depth  set depth of selected items to depth\n";
  print STDERR "  -trans=dx,dy  translate selected nodes by dx/dy\n";
  print STDERR "  -delete       delete selected items\n";
  print STDERR "                                    \n";
  print STDERR "  -invert       invert selection - modify outer area\n";
  print STDERR "  -exclude      exclude border line items from selection\n";
  print STDERR "  -preserve     preserves order of items\n";
  print STDERR "  -depth_invert inverts depth values (neg to pos etc.)\n";
  print STDERR "  -unset_depth  deletes depth values\n";
}

sub Usage {
  print STDERR "Usage: grd_modify.pl [-h|-help] [-options] grid [line]\n";
}

#----------------------------------------------------------

sub loop_on_nodes {

  my ($grid,$flag,$proc) = @_;

  $grid->make_used();
  my $nodes = $grid->get_nodes();

  foreach my $node (values %$nodes) {
    my $number = $node->{number};
    my $f = $$flag{$number};
    if( $f and $::inside or not $f and not $::inside ) {
      &$proc($grid,$node);
    }
  }
}

sub loop_on_elements {

  my ($grid,$flag,$proc) = @_;

  my $elems = $grid->get_elems();

  loop_on_items($grid,$flag,$proc,$elems);
}
 
sub loop_on_lines {

  my ($grid,$flag,$proc) = @_;

  my $lines = $grid->get_lines();

  loop_on_items($grid,$flag,$proc,$lines);
}

sub loop_on_items {

  my ($grid,$flag,$proc,$items) = @_;

  my $elems = $grid->get_elems();
  foreach my $item (values %$items) {
    my $vert = $item->{vert};
    my $nvert = $item->{nvert};
    my $number = $item->{number};
    my $n = 0;
    foreach my $node (@$vert) {
      $n++ if $$flag{$node};
    }

#    if( $n ) {
#	print "flagged: element $number ($n)\n";
#    }

    if( $::inside ) {
      if( $::exclude ) {
        if( $n == $nvert ) {
	  &$proc($grid,$item);
	}
      } else {
        if( $n ) {
	  &$proc($grid,$item);
	}
      }
    } else {
      if( $::exclude ) {
        if( $n == 0 ) {
	  &$proc($grid,$item);
	}
      } else {
        if( $n < $nvert ) {
	  &$proc($grid,$item);
	}
      }
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
  } elsif( $::type != $::flag ) {
    $elem->{type} = $::type;
  } elsif( $::delete ) {
    $grid->delete_elem($elem);
  } elsif( $::depth != $::flag ) {
    $elem->{h} = $::depth;
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
  } elsif( $::type != $::flag ) {
    $line->{type} = $::type;
  } elsif( $::delete ) {
    $grid->delete_line($line);
  } elsif( $::depth != $::flag ) {
    $line->{h} = $::depth;
  } elsif( $::depth_invert ) {
    $line->{h} = -$line->{h};
  } elsif( $::unset_depth ) {
    $line->{h} = $::depth_flag;
  }
}

#-----------------------------------------------------------------

sub print_flaged_nodes {

  my ($flag) = @_;

  foreach my $number (keys %$flag) {
    print "$number\n";
  }
}

sub flag_nodes {

  my ($grid,$grline,$flag) = @_;

  my $nodes = $grid->get_nodes();
  foreach my $node (values %$nodes) {
    my $number = $node->{number};
    my $x = $node->{x};
    my $y = $node->{y};
    if( $grline->in_line($x,$y) ) {
      $$flag{$number} = 1;
    }
  }

  #return $flag;
}

#-----------------------------------------------------------------

sub make_grdline {

  my ($gline,$line) = @_;

  my ($x,$y) = make_xy_array($gline,$line);
  my $grline = new grdline;
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

