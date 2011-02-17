#!/usr/bin/perl -w
#
# grd utility routines
#
##############################################################
#
# version 1.7
#
# 19.08.2005		unify_nodes, if defined $depth
# 24.08.2005		connect_lines, split_line, contains_node
# 26.08.2005		make_xy, node_info, elem_info
# 20.10.2005		version control introduced
# 01.03.2010		preserve order of written items if needed
# 07.10.2010		new routine get_xy_minmax()
# 10.02.2011		new routine make_central_point()
#
##############################################################
#
# example of usage:
#
# use lib "$ENV{HOME}/lib/perl";
# use grd;
# 
# my $grid = new grd;
# my $infile = $ARGV[0];
# 
# $grid->readgrd($infile);
# $grid->writegrd("out.grd");
# 
# to iterate over say the list of elements:
#
# my $elems = $grid->get_elems();
# foreach my $item (values %$elems) {
#   my $number = $item->{number};
#   ...
# }
#
# or
#
# foreach my $item ( $grid->get_elem_list() ) {
#   my $number = $item->{number};
#   ...
# }
#
##############################################################

use strict;

package grd;

##############################################################

sub new
{
    my $self;

    my %nodes = ();
    my %elems = ();
    my %lines = ();
    my @comms = ();

    $self =	{
	    		 file		=>	undef
	    		,outfile	=>	undef
	    		,outhandle	=>	*STDOUT
			,nodes		=>	\%nodes
			,elems		=>	\%elems
			,lines		=>	\%lines
			,comms		=>	\@comms
			,nodes_ordered	=>	[]
			,elems_ordered	=>	[]
			,lines_ordered	=>	[]
			,nnmax		=>	0
			,nemax		=>	0
			,nlmax		=>	0
			,verbose	=>	1
			,preserve_order	=>	0
		};

    bless $self;
    return $self;
}

###############################################################################

sub set_outfile {

  my ($self,$outfile) = @_;

  if( $self->{outhandle} and $self->{outhandle} ne *STDOUT ) {
      close($self->{outhandle});
  }

  my $handle;
  if( $outfile ) {
    my $file = make_grdfile($outfile);
    open($handle,">$outfile") || die "Cannot open outfile: $file\n";
  } else {
    $handle = *STDOUT;
  }

  $self->{outfile} = $outfile;
  $self->{outhandle} = $handle;
}

sub writegrd {

  my ($self,$file) = @_;

  $self->set_outfile($file);
  my $fh = $self->{outhandle};
  $file = $self->{outfile};

  if( $self->{verbose} ) {
    print STDERR "writing file: $file\n" if $file;
  }

  my ($items,@keys);

  print $fh "\n";

  foreach my $line (@{$self->{comms}}) {
    print $fh "$line\n";
  }
  print $fh "\n";
    
  $items = $self->{nodes};
  @keys = keys %$items;
  @keys = sort {$a<=>$b} @keys;
  @keys = @{$self->{nodes_ordered}} if $self->{preserve_order};
  foreach my $number (@keys) {
      $self->write_node($items->{$number}) if $items->{$number};
  } 
  print $fh "\n";

  $items = $self->{elems};
  @keys = keys %$items;
  @keys = sort {$a<=>$b} @keys;
  @keys = @{$self->{elems_ordered}} if $self->{preserve_order};
  foreach my $number (@keys) {
      $self->write_elem($items->{$number});
  } 
  print $fh "\n";

  $items = $self->{lines};
  @keys = keys %$items;
  @keys = sort {$a<=>$b} @keys;
  @keys = @{$self->{lines_ordered}} if $self->{preserve_order};
  foreach my $number (@keys) {
      $self->write_line($items->{$number});
  } 
  print $fh "\n";

}

sub write_nodes {		# helper routine

  my ($fh,$ra) = @_;
  my $j = 0;

  foreach (@$ra) {
    $j++;
    print $fh " $_";
    print $fh "\n" if $j%10 == 0;
  }
}

sub print_node
{
    my ($self,$number) = @_;

    my $node = $self->{nodes}->{$number};

    $self->write_node($node);
}

sub write_node
{
    my ($self,$item) = @_;

    my $fh = $self->{outhandle};

    print $fh "1 $item->{number} $item->{type}";
    print $fh " $item->{x} $item->{y}";
    print $fh " $item->{h}" if defined $item->{h};
    print $fh "\n";
}

sub write_elem
{
    my ($self,$item) = @_;

    my $fh = $self->{outhandle};
    my $nvert = $item->{nvert};

    print $fh "2 $item->{number} $item->{type}";
    print $fh " $item->{nvert}";
    print $fh "\n" if $nvert > 3;
    &write_nodes($fh,$item->{vert});
    print $fh "\n" if $nvert > 3 and $item->{h};
    print $fh " $item->{h}" if defined $item->{h};
    print $fh "\n";
} 

sub write_line
{
    my ($self,$item) = @_;

    my $fh = $self->{outhandle};

    print $fh "3 $item->{number} $item->{type} $item->{nvert}\n";
    &write_nodes($fh,$item->{vert});
    print $fh " $item->{h}" if defined $item->{h};
    print $fh "\n\n";
}

###############################################################################

sub readgrd {

  my ($self,$file) = @_;

  return unless $file;

  $file = make_grdfile($file);
  open(FILE,"$file") || die "Cannot open file: $file\n";

  $self->{file} = $file;
  if( $self->{verbose} ) {
    print STDERR "reading file: $file\n";;
  }

  while( $_ = $self->nextitem ) {

	#print "readgrd: $_\n";

        next if /^\s*$/;
	my @f = split;

	my $item = $f[0];

	if( $item == 0 ) {
		$self->insert_comment($_);
	} elsif( $item == 1 ) {
		$self->insert_node(\@f);
	} elsif( $item == 2 ) {
		$self->insert_elem(\@f);
	} elsif( $item == 3 ) {
		$self->insert_line(\@f);
	} else {
		die "Unknown item: $_\n";
	}
  }

  close(FILE);
}

sub make_grdfile {

  my $name = shift;
  my $file = "";

  if( $name =~ /\.grd\s*$/i ) {
    $file = $name;
  } else {
    $file = "$name.grd";
  }

  return $file;
}

###################################

sub nextitem {

  my ($self) = @_;

  my $line = &getline;
  return undef unless defined($line);

  #print "nextitem (line): $line";

  while( my $newline = &getline ) {
    #print "nextitem (newline): $newline";
    if( $newline =~ /^\s*$/ ) {	#empty line
	last;
    } elsif( $newline =~ /^\s+/ ) {		#conti line
	$line .= " $newline";
    } else {
	&ungetline($newline);		#save for next call
	last;
    }
  }

  $line =~ s/\n/ /g;

  return $line;
}

###################################

sub getline {

  if( @grd::oldlines ) {
    return shift(@grd::oldlines);
  } else {
    return <FILE>;
  }
}

sub ungetline {

  push(@grd::oldlines,$_[0]);
}

###################################

sub insert_comment 
{
    my ($self,$line) = @_;

    my $rcomms = $self->{comms};
    push(@$rcomms,$line);
}

sub insert_node 
{
    my ($self,$ritems) = @_;

    my %item = ();

    my $number = $$ritems[1];
    my $rnodes = $self->{nodes};

    $item{number} = $number;
    $item{type}   = $$ritems[2];
    $item{x}      = $$ritems[3];
    $item{y}      = $$ritems[4];
    $item{h}      = $$ritems[5];

    $item{used}   = 0;
    $self->{nnmax} = $number if $number > $self->{nnmax};

    $rnodes->{$number} = \%item;
    push(@{$self->{nodes_ordered}},$number);
}

sub insert_elem 
{
    my ($self,$ritems) = @_;

    my %item = ();

    my $number = $$ritems[1];
    my $nvert  = $$ritems[3];
    my $relems = $self->{elems};

    $item{number} = $number;
    $item{type}   = $$ritems[2];
    $item{nvert}  = $nvert;
    $item{vert}   = &get_vertices($ritems);
    $item{h}      = $$ritems[4+$nvert];

    $self->{nemax} = $number if $number > $self->{nemax};

    $relems->{$number} = \%item;
    push(@{$self->{elems_ordered}},$number);
}

sub insert_line 
{
    my ($self,$ritems) = @_;

    my %item = ();

    my $number = $$ritems[1];
    my $nvert  = $$ritems[3];
    my $rlines = $self->{lines};

    $item{number} = $number;
    $item{type}   = $$ritems[2];
    $item{nvert}  = $nvert;
    $item{vert}   = &get_vertices($ritems);
    $item{h}      = $$ritems[4+$nvert];

    $self->{nlmax} = $number if $number > $self->{nlmax};

    $rlines->{$number} = \%item;
    push(@{$self->{lines_ordered}},$number);
}

###############################################################################

sub node_info
{
    my ($self,$item,$text) = @_;

    print STDERR "$text\n";
    if( $item ) {
      print STDERR "$item->{number} $item->{type} $item->{x} $item->{y}\n";
    } else {
      print STDERR "no such item\n";
    }
}

sub elem_info
{
    my ($self,$item,$text) = @_;

    print STDERR "$text\n";
    if( $item ) {
      print STDERR "$item->{number} $item->{type} $item->{nvert}\n";
    } else {
      print STDERR "no such item\n";
    }
}

###############################################################################

sub by_number {
  $a<=>$b;
}

sub compress_node_numbers
{
    my ($self) = @_;

    my $nodes = $self->{nodes};
    my @numbers = sort by_number keys(%$nodes);
    my %intern = ();
    my %new_nodes = ();
    my $items;

    my $i = 1;
    foreach my $number (@numbers) {
      $intern{$number} = $i;
      my $node = $nodes->{$number};
      $node->{number} = $i;
      $new_nodes{$i} = $node;
      $i++;
    }
    $self->{nodes} = \%new_nodes;

    $items = $self->{elems};
    &compress_node_list($items,\%intern);

    $items = $self->{lines};
    &compress_node_list($items,\%intern);
}

sub compress_node_list
{

    my ($items,$intern) = @_;

    my @keys = keys %$items;
    foreach my $number (@keys) {
      my $item = $items->{$number};
      my $node_list = $item->{vert};
      my @new_nodes = ();
      foreach my $node (@$node_list) {
	push(@new_nodes,$$intern{$node});
      }
      $item->{vert} = \@new_nodes;
    }
}

sub get_xy_minmax
{
    my $self = shift;

    my $nodes = $self->{nodes};

    my $xmin = 1.e+30;
    my $ymin = 1.e+30;
    my $xmax = -1.e+30;
    my $ymax = -1.e+30;

    foreach my $item (values %$nodes) {
      my $x = $item->{x};
      my $y = $item->{y};
      $xmin = $x if $x < $xmin;
      $ymin = $y if $y < $ymin;
      $xmax = $x if $x > $xmax;
      $ymax = $y if $y > $ymax;
   }

  return ($xmin,$ymin,$xmax,$ymax);
}

sub make_central_point 
{
    my ($self,$item) = @_;	#item must be element or line

    my $xm = 0;
    my $ym = 0;

    my $verts = $item->{vert};
    foreach my $nnode (@$verts) {
      my $node = $self->get_node($nnode);
      my $x = $node->{x};
      my $y = $node->{y};
      $xm += $x;
      $ym += $y;
    }
    my $nverts = $item->{nvert};
    if( $verts->[0] == $verts->[-1] ) {
      $nverts -= 1;
    }

    $xm /= $nverts;
    $ym /= $nverts;

    return ($xm,$ym);
}


###############################################################################

sub get_node_list { return $_[0]->get_item_list("nodes"); }
sub get_elem_list { return $_[0]->get_item_list("elems"); }
sub get_line_list { return $_[0]->get_item_list("lines"); }

sub get_item_list
{
    my ($self,$type) = @_;

    my $items = $self->{$type};
    return values %$items
}

#----------

sub get_nodes { return $_[0]->get_items("nodes"); }
sub get_elems { return $_[0]->get_items("elems"); }
sub get_lines { return $_[0]->get_items("lines"); }

sub get_items
{
    my ($self,$type) = @_;

    return $self->{$type};
}

#----------

sub exists_node { return $_[0]->exists_item($_[1],"nodes"); }
sub exists_elem { return $_[0]->exists_item($_[1],"elems"); }
sub exists_line { return $_[0]->exists_item($_[1],"lines"); }

sub exists_item
{
    my ($self,$number,$type) = @_;

    my $items = $self->{$type};
    if( exists $items->{$number} ) {
      return 1;
    } else {
      return 0;
    }
}

#----------

sub get_node { return $_[0]->get_item($_[1],"nodes"); }
sub get_elem { return $_[0]->get_item($_[1],"elems"); }
sub get_line { return $_[0]->get_item($_[1],"lines"); }

sub get_item
{
    my ($self,$number,$type) = @_;

    my $items = $self->{$type};
    unless( exists $items->{$number} ) {
      die "get_item: Cannot find number $number of type $type\n";
    }
    return $items->{$number};
}

#----------
#
# $grid->delete_node($node)    where $node is ref to node item
# $grid->delete_node(511)      directly with number

sub delete_node { return $_[0]->delete_item($_[1],"nodes"); }
sub delete_elem { return $_[0]->delete_item($_[1],"elems"); }
sub delete_line { return $_[0]->delete_item($_[1],"lines"); }

sub delete_item
{
    my ($self,$item,$type) = @_;

    my $items = $self->{$type};
    $item = $self->get_item($item,$type) unless ref($item); #if number, get item
    my $number = $item->{number};

    delete $$items{$number};
}

#----------

sub clone_node { return $_[0]->clone_item($_[1],"nodes"); }
sub clone_elem { return $_[0]->clone_item($_[1],"elems"); }
sub clone_line { return $_[0]->clone_item($_[1],"lines"); }

sub clone_item
{
    my ($self,$item,$type) = @_;

    my $items = $self->{$type};
    my $number = $item->{number};
    if( exists $items->{$number} ) {
      print STDERR "clone_item: Item $number of type $type exists already\n";
    } else {
      $items->{$number} = $item;
    }
}

###############################################################################

sub unify_nodes {	#unifies nodes -> node n2 is deleted

  my ($self,$n1,$n2) = @_;

  $n1 = $n1->{number} if ref($n1);	# get number if item
  $n2 = $n2->{number} if ref($n2);	# get number if item

  $self->delete_node($n2);
  $self->substitute_node($n1,$n2);
}

sub substitute_node {

  my ($self,$n1,$n2) = @_;

  # substitute n2 with n1

  $self->substitute_vert($self->{elems},$n1,$n2);
  $self->substitute_vert($self->{lines},$n1,$n2);
}

sub substitute_vert {

  my ($self,$items,$n1,$n2) = @_;

  # substitute n2 with n1

  foreach my $item (values %$items) {
    my $vert = $item->{vert};
    foreach my $n (@$vert) {
	$n = $n1 if $n == $n2;
    }
  }
}

###################################

sub make_node
{
    my ($self,$number,$type,$x,$y,$depth) = @_;

    my %item = ();

    if( $number <= 0 ) {
	$self->{nnmax}++;
	$number = $self->{nnmax};
    }

    $item{number} = $number;
    $item{type}   = $type;
    $item{x}      = $x;
    $item{y}      = $y;
    $item{h}      = $depth if defined $depth;

    $self->{nodes}->{$number} = \%item;

    return \%item;
}

sub make_elem
{
    my ($self,$number,$type,$depth,@vert) = @_;

    my %item = ();

    if( $number <= 0 ) {
	$self->{nemax}++;
	$number = $self->{nemax};
    }

    #print STDERR "make_elem: $vert[0] $vert[1] $vert[2]\n";

    $item{number} = $number;
    $item{type}   = $type;
    $item{nvert}  = @vert;
    $item{vert}   = \@vert;
    $item{h}      = $depth if defined $depth;

    $self->{elems}->{$number} = \%item;

    return \%item;
}

sub make_line
{
    my ($self,$number,$type,$depth,@vert) = @_;

    my %item = ();

    if( $number <= 0 ) {
	$self->{nlmax}++;
	$number = $self->{nlmax};
    }

    $item{number} = $number;
    $item{type}   = $type;
    $item{nvert}  = @vert;
    $item{vert}   = \@vert;
    $item{h}      = $depth if defined $depth;

    $self->{lines}->{$number} = \%item;

    return \%item;
}

###################################

sub make_xy
{
	my ($self,$item) = @_;

	my @x = ();
	my @y = ();
	my $v = $item->{vert};

	foreach my $n (@$v) {
	  my $node = $self->get_node($n);
	  push(@x,$node->{x});
	  push(@y,$node->{y});
	}

	return (\@x,\@y);
}

sub contains_node
{
  my ($self,$item,$n) = @_;

  my $v = $item->{vert};

  foreach my $node (@$v) {
    return 1 if $node == $n;
  }
  return 0;
}

sub is_closed
{
  my ($self,$line) = @_;

  my $vert = $line->{vert};

  if( $$vert[0] == $$vert[-1] ) {
	return 1;
  } else {
	return 0;
  }
}

sub close_line
{
  my ($self,$line) = @_;

  my $vert = $line->{vert};

  if( $$vert[0] != $$vert[-1] ) {
    push(@$vert,$$vert[0]);
    $line->{nvert}++;
  }
}

sub connect_lines
{
  # we try to use node, but if node is wrong we connect anyway

  my ($self,$line1,$line2,$node) = @_;

  my $vert1 = $line1->{vert};
  my $vert2 = $line2->{vert};
  my $n = $node->{number} if $node;

  my @new = ();
  my @v1 = ();
  my @v2 = ();

  if( $$vert1[-1] == $$vert2[0] ) {
    @v1 = @$vert1;
    @v2 = @$vert2;
  } elsif( $$vert1[0] == $$vert2[0] ) {
    @v1 = reverse @$vert1;
    @v2 = @$vert2;
  } elsif( $$vert1[-1] == $$vert2[-1] ) {
    @v1 = @$vert1;
    @v2 = reverse @$vert2;
  } elsif( $$vert1[0] == $$vert2[-1] ) {
    @v1 = reverse @$vert1;
    @v2 = reverse @$vert2;
  } else {
    print STDERR  "*** Cannot connect lines $line1 and $line2: no node in common\n";
    return;
  }

  if( $node and $v2[0] != $n ) {	#node is given and link not ok
    if( $v1[0] == $n and $v2[-1] == $n ) {
	my @aux = @v1;
	@v1 = @v2;
	@v2 = @aux;
    }
  }
  push(@new,@v1);
  shift(@v2);
  push(@new,@v2);

  $line1->{vert} = \@new;
  $line1->{nvert} = @new;

  $self->delete_line($line2);

  return $line1;
}

sub split_line
{
  my ($self,$line,$node) = @_;

  my $number = $node->{number};
  my $vert = $line->{vert};
  my @v1 = ();
  my @v2 = ();
  my $act = \@v1;

  foreach my $n (@$vert) {
	push(@$act,$n);
	if( $n == $number ) {
	  $act = \@v2;
	  push(@$act,$n);
	}
  }
	  
  if( scalar(@v1) < 2 or scalar(@v2) < 2 ) {
    print STDERR "*** cannot split line on this point\n";
    return ($line,"");
  }

  $line->{nvert} = scalar(@v1);
  $line->{vert} = \@v1;

  my $line_new = $self->make_line(0,$line->{type},$line->{depth},@v2);

  return ($line,$line_new);
}

###################################

sub get_vertices
{
    my $ritems = shift;

    my @array = @$ritems;
    my @new = ();

    shift(@array);
    shift(@array);
    shift(@array);
    my $nvert = shift(@array);

    foreach my $vert (@array) {
	push(@new,$vert);
	$nvert--;
	last unless $nvert;
    }

    die "*** Not enough vertices... @$ritems\n" if $nvert;

    return \@new;
}

###################################

sub delete_unused
{
    my ($self) = @_;

    $self->make_used();

    my $nodes = $self->{nodes};
    foreach my $node (values %$nodes) {
      unless( $node->{used} ) {
        $self->delete_node($node);
      }
    }
}

sub make_used
{
    my ($self) = @_;

    my $nodes = $self->{nodes};
    foreach my $node (values %$nodes) {
      $node->{used} = 0;
    }

    my $elems = $self->{elems};
    foreach my $item (values %$elems) {
      $self->inc_used($item->{vert});
    }

    my $lines = $self->{lines};
    foreach my $item (values %$lines) {
      $self->inc_used($item->{vert});
    }
}

sub inc_used
{
    my ($self,$verts) = @_;
    my $nodes = $self->{nodes};

    foreach my $number (@$verts) {
      my $node = $nodes->{$number};
      #print STDERR 
      $node->{used}++;
    }
}

###################################

sub set_verbose
{
    my ($self,$verbose) = @_;

    $self->{verbose} = $verbose;
}

sub set_preserve_order
{
    my ($self,$preserve_order) = @_;

    $self->{preserve_order} = $preserve_order;
}

###################################

sub test_grd
{
    my @files = @_;

    my $grid = new grd;

    print "reading...\n";
    $grid->readgrd($files[0]);
    print "writing...\n";
    $grid->writegrd;
}

###################################
#&test_grd(@ARGV);
###################################
1;
###################################

