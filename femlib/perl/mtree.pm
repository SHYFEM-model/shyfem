#!/usr/bin/perl -w
#
##############################################################
#
# version 1.0
#
# 01.12.2011            written from scratch
#
##############################################################
#
# implements multiple tree with contained feature
#
#	objects contain other objects -> a hirarchy is constructed
#	two functions have to be provided:
#	inside(item1,item2)	is true if item2 is inside item1
#	get_id(item)		returns integer id of item
#
# usage:
#
# use mtree;
#
# my $tree = new mtree(\&inside,\&get_id);
#
# or
#
# my $tree = new mtree;
# $tree->set_functions(\&inside,\&get_id);
#
# while(...) {
#   $tree->insert($item)
# }
# $tree->print()
#
##############################################################
 
use strict;
 
package mtree;

sub new
{
    my ($self,$inside,$get_id) = @_;
    #my $self;
 
    $self =     {
                         root        =>      {}
                        ,inside      =>      $inside
                        ,get_id      =>      $get_id
                };
 
    bless $self;
    return $self;
}
 
sub set_functions {

  my ($self,$inside,$get_id) = @_;

  $self->{inside} = $inside;
  $self->{get_id} = $get_id;
}

##############################################################

sub insert {

  my ($self,$item) = @_;

  $self->insert_item($self->{root},$item);
}

sub insert_item {

  my ($self,$node,$item) = @_;

  my $inside = $self->{inside};

  my $children = $node->{children};

  my $siblings = [];
  my $childs = [];

  if( $children ) {
    foreach my $child (@$children) {
      my $citem = $child->{info};
      if( $inside->($item,$citem) ) {		# child inside new item
	push(@$childs,$child);
      } elsif( $inside->($citem,$item) ) {	# new item inside child
        $self->insert_item($child,$item);
	return;
      } else {						# unrelated
	push(@$siblings,$child);
      }
    }
  }

  my $new_node = $self->new_node($node,$childs,$item);	#parent, children, info

  foreach my $child (@$childs) {
    $child->{parent} = $new_node;
  }

  push(@$siblings,$new_node);
  $node->{children} = $siblings;
}

sub new_node {

  my ($self,$parent,$children,$info) = @_;

  my $new_node = {};
  $new_node->{parent} = $parent;
  $new_node->{children} = $children;
  $new_node->{info} = $info;

  return $new_node;
}

#---------------------------------------------------

sub get_child_items {

  my ($self,$children) = @_;

  my @list;

  foreach my $child (@$children) {
    push(@list,$child->{info});
  }

  return \@list;
}

#---------------------------------------------------

sub apply {

  my ($self,$func) = @_;

  $self->apply_node($self->{root},$func);
}

sub apply_node {

  my ($self,$node,$func) = @_;

  my $children = $node->{children};
  if( $children ) {
    foreach my $child (@$children) {
      $self->apply_node($child,$func);
    }
  }

  my $inode = $node->{info};
  return unless $inode;

  my $parent = $node->{parent};
  my $iparent = $parent->{info} if $parent;

  my $ichilds = $self->get_child_items($children);

  $func->($inode,$iparent,$ichilds);
}

#---------------------------------------------------

sub print {

  my ($self) = @_;

  $self->print_node($self->{root});
}

sub print_node {

  my ($self,$node) = @_;

  my $space = $::space;
  $::space .= "  ";

  my $line = $self->get_node_info($node);
  print "$::space $line\n";

  my $children = $node->{children};
  if( $children ) {
    foreach my $child (@$children) {
      $self->print_node($child);
    }
  }
    
  $::space = $space;
}

sub get_node_info {

  my ($self,$node) = @_;

  my $get_id = $self->{get_id};

  my $info = $node->{info};
  my $n = $get_id->($info) if $info;
  $n = 0 unless $n;

  my $parent = $node->{parent};
  my $pinfo = $parent->{info} if $parent;;
  my $pn = $get_id->($pinfo) if $pinfo;
  $pn = 0 unless $pn;

  my $cline = "";
  my $children = $node->{children};
  $children = [] unless $children;
  foreach my $child (@$children) {
    my $cinfo = $child->{info};
    #my $cn = $cinfo->{number};
    my $cn = $get_id->($cinfo);
    $cline .= " $cn";
  }

  return "line $n (parent: $pn; children: $cline)"
}
  
##############################################################
1;
##############################################################

