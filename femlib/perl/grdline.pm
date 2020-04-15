#!/usr/bin/perl -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
##############################################################
#
# version 1.4.0
#
# 26.08.2005            print_info, bug in make_unique
# 20.10.2005            version control introduced
# 01.12.2011            new functionality
# 18.02.2014            bug fix in is_convex()
#
##############################################################
#
# usage:
#
# use grdline;
#
# my $grdl = new grdline;
# $grdl->set_line($x,$y);
#
# if package grd.pm is used you can do:
#
# ($x,$y) = $grid->make_xy($litem);
# $grdl->set_line($x,$y);
#
##############################################################
 
use strict;
 
package grdline;
 
##############################################################
 
sub new
{
    my $self;
 
    my %x = ();
    my %y = ();
 
    $self =     {
                         n           =>      0
                        ,x           =>      \%x
                        ,y           =>      \%y
                        ,closed      =>      -1
                        ,convex      =>      -1
                        ,xmin        =>      -1
                        ,xmax        =>      -1
                        ,ymin        =>      -1
                        ,ymax        =>      -1
                        ,area        =>      -1
                };
 
    bless $self;
    return $self;
}
 
##############################################################

sub reset_line {

  my ($self) = @_;

  $self->{closed} = -1;
  $self->{convex} = -1;
  $self->{area} = -1;
}

sub set_line {

  my ($self,$x,$y) = @_;

  $self->reset_line();

  my $nx = @$x;
  my $ny = @$y;

  die "Coordinates of different length: $nx $ny\n" if $nx != $ny;

  $self->{n} = $nx;

  my @xnew = @$x;
  my @ynew = @$y;

  $self->{x} = \@xnew;
  $self->{y} = \@ynew;

  $self->make_unique();
  $self->is_closed();
  $self->set_area();
  $self->is_convex();
  $self->set_xy_min_max();
}

sub print_info {

  my ($self,$text,$verbose) = @_;

  print STDERR "$text\n";
  print STDERR "n: $self->{n}   closed: $self->{closed}" .
    "   convex: $self->{convex}\n";
  print STDERR "xy-min/max: $self->{xmin} $self->{ymin}" .
    " $self->{xmax} $self->{ymax}\n";
}

sub is_convex {

  my ($self) = @_;

  return $self->{convex} if $self->{convex} != -1;

  my $area = $self->{area};
  $self->{convex} = 0;

  my $n = $self->{n} - 1;
  my $x = $self->{x};
  my $y = $self->{y};

  my ($xl,$yl);
  my ($xm,$ym) = ($$x[$n-1],$$y[$n-1]);
  my ($xn,$yn) = ($$x[$n],$$y[$n]);

  for (my $i=0;$i<=$n;$i++) {
    ($xl,$yl) = ($xm,$ym);
    ($xm,$ym) = ($xn,$yn);
    ($xn,$yn) = ($$x[$i],$$y[$i]);

    if( $area > 0 ) {
      return 0 if not lefton($xl,$yl,$xm,$ym,$xn,$yn);
    } elsif( $area < 0 ) {
      return 0 if not righton($xl,$yl,$xm,$ym,$xn,$yn);
    } else {
      return 0;
    }
  }

  $self->{convex} = 1;
  return 1;
}

sub in_convex {

  my ($self,$x0,$y0) = @_;

  my $n = $self->{n} - 1;
  my $x = $self->{x};
  my $y = $self->{y};

  my ($xm,$ym);
  my ($xn,$yn) = ($$x[$n],$$y[$n]);

  for (my $i=0;$i<=$n;$i++) {
    ($xm,$ym) = ($xn,$yn);
    ($xn,$yn) = ($$x[$i],$$y[$i]);

    #print STDERR "in_convex: $n $i $xm,$ym,$xn,$yn,$x0,$y0\n";
    return 0 unless( lefton($xm,$ym,$xn,$yn,$x0,$y0) );
  }

  return 1;
}

sub make_closed {
  my ($self) = @_;
  $self->{closed} = 1;
}

sub is_closed {

  my ($self) = @_;

  return $self->{closed} if $self->{closed} != -1;

  my $n = $self->{n} - 1;
  my $x = $self->{x};
  my $y = $self->{y};

  if( $$x[0] != $$x[$n] or $$y[0] != $$y[$n] ) {
    $self->{closed} = 0;
    return 0;
  }

  $self->{n} = $n;
  pop(@$x);
  pop(@$y);

  $self->{closed} = 1;
  return 1;
}

sub make_unique {

  my ($self) = @_;

  my $n = $self->{n};
  my $x = $self->{x};
  my $y = $self->{y};

  my @newx = ($x->[0]);
  my @newy = ($y->[0]);

  for( my $i=1; $i<$n; $i++) {
    if( $$x[$i] != $$x[$i-1] or $$y[$i] != $$y[$i-1] ) {
      push(@newx,$$x[$i]);
      push(@newy,$$y[$i]);
    } else {
      print STDERR "*** eliminating double point...\n";
    }
  }

  $self->{n} = @newx;
  $self->{x} = \@newx;
  $self->{y} = \@newy;
}

sub set_area {
 
  my $self = shift;

  $self->{area} = $self->area();
}

sub area {
 
  my $self = shift;

  my $n = $self->{n}-1;
  my $x = $self->{x};
  my $y = $self->{y};

  my ($xc,$yc) = $self->get_center_point();
  my $area = 0.;
 
  my ($xm,$ym);
  my ($xn,$yn) = ($$x[$n],$$y[$n]);

  for (my $i=0;$i<=$n;$i++) {
    ($xm,$ym) = ($xn,$yn);
    ($xn,$yn) = ($$x[$i],$$y[$i]);

    $area += areat($xm,$ym,$xn,$yn,$xc,$yc);
  }

  return $area;
}

sub get_center_point {
 
  my $self = shift;
  my $n = $self->{n};
  my $x = $self->{x};
  my $y = $self->{y};

  my ($xc,$yc) = (0,0);
 
  for( my $i=0;$i<$n;$i++ ) {
    $xc += $$x[$i];
    $yc += $$y[$i];
  }
 
  return ($xc/$n,$yc/$n);
}

#--------------------------------------------------------------------

sub areat {

  my ($x1,$y1,$x2,$y2,$x3,$y3) = @_;

  return 0.5 * ( ($x2-$x1) * ($y3-$y1) - ($x3-$x1) * ($y2-$y1) );
}

sub lefton  { return ( areat(@_) >= 0 ); }
sub righton { return ( areat(@_) <= 0 ); }

sub angle {

  my ($x0,$y0,$xa,$ya,$xb,$yb) = @_;

  # angle between a-0-b (angle is on 0)

  my $dxa = $xa - $x0;
  my $dya = $ya - $y0;
  my $dxb = $xb - $x0;
  my $dyb = $yb - $y0;

  my $moda = sqrt( $dxa*$dxa + $dya*$dya );
  my $modb = sqrt( $dxb*$dxb + $dyb*$dyb );
  my $mod = $moda * $modb;

  if( $mod <= 0. ) {
    return 0.;
  }

  my $cos = (  $dxa * $dxb + $dya * $dyb ) / $mod;
  my $sin = ( -$dya * $dxb + $dxa * $dyb ) / $mod;

  if( $cos > 1. ) {
    $cos = 1.;
  }

  my $alpha = atan2( sqrt(1-$cos*$cos) , $cos );
  $alpha = -$alpha if( $sin < 0. );

  return $alpha;
}

sub get_min_max {
 
  my $a = shift;
 
  my ($min,$max);
 
  $min = $max = $$a[0];
  foreach my $c (@$a) {
    $max = $c if $c > $max;
    $min = $c if $c < $min;
  }
 
  return ($min,$max);
}

#--------------------------------------------------------------------

sub set_xy_min_max {
 
  my ($self) = @_;
 
  ($self->{xmin},$self->{xmax}) = get_min_max($self->{x});
  ($self->{ymin},$self->{ymax}) = get_min_max($self->{y});
}

sub in_min_max {

  my ($self,$x0,$y0) = @_;

  if( $self->{xmin} <= $x0 and $x0 <= $self->{xmax} ) {
    if( $self->{ymin} <= $y0 and $y0 <= $self->{ymax} ) {
	return 1;
    }
  }
  return 0;
}

sub line_in_line {

  my ($self,$grdl) = @_;

  my $x = $grdl->{x};
  my $y = $grdl->{y};

  my $x0 = $x->[0];
  my $y0 = $x->[0];

  return $self->in_line($x0,$y0);
}

sub in_line {

  my ($self,$x0,$y0) = @_;

  #print STDERR "in_line: $self->{n} $self->{convex} $self->{closed}\n";

  if( $self->{convex} ) {
    return $self->in_convex($x0,$y0);
  }

  return 0 unless $self->in_min_max($x0,$y0);

  return $self->in_any_line($x0,$y0);
}

sub in_any_line {

  my ($self,$x0,$y0) = @_;

  my $n = $self->{n} - 1;
  my $x = $self->{x};
  my $y = $self->{y};

  my $angle;
  my $pi = 4*atan2(1,1);

  my ($xm,$ym);
  my ($xn,$yn) = ($$x[$n],$$y[$n]);

  for (my $i=0;$i<=$n;$i++) {
    ($xm,$ym) = ($xn,$yn);
    ($xn,$yn) = ($$x[$i],$$y[$i]);

    $angle += angle($x0,$y0,$xm,$ym,$xn,$yn);
  }

  if( abs($angle) < $pi ) {
    return 0;
  } else {
    return 1;
  }
}

#--------------------------------------------------------------------

###################################
1;
###################################

