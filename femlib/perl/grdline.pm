#!/usr/bin/perl -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
##############################################################
#
# version 1.6.0
#
# 26.08.2005            print_info, bug in make_unique
# 20.10.2005            version control introduced
# 01.12.2011            new functionality
# 18.02.2014            bug fix in is_convex()
# 05.07.2021            new functionality
# 13.07.2021            bug fix, documentation
# 15.06.2022            new function invert_line() (must be set individually)
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
# ($xr,$yr) = $grid->make_xy($litem);
# $grdl->set_line($xr,$yr);		# $xr,$yr are references to array
#
# implemented calls:
#
# $grdl->print_info($text);		# $text is explicative string
# $grdl->make_closed();
# $area = $grdl->get_area();
# ($xc,$yc) = $grdl->get_center_point();
# ($xmin,$ymin,$xmax,$ymax) = $grdl->get_xy_min_max();
# 
# $flag = $grdl->is_convex();
# $flag = $grdl->is_closed();
# $flag = $grdl->in_line($x,$y);	# $x,$y is point to be checked
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
                        ,setup       =>      0
                        ,closed      =>      -1
                        ,convex      =>      -1
                        ,latlon      =>      -1
                        ,must_invert =>       0
                        ,xmin        =>      -1
                        ,xmax        =>      -1
                        ,ymin        =>      -1
                        ,ymax        =>      -1
                        ,xc          =>      -1
                        ,yc          =>      -1
                        ,area        =>      -1
                };
 
    bless $self;
    return $self;
}
 
##############################################################

sub reset_line {

  my ($self) = @_;

  $self->{setup} = 0;
}

sub set_line {

  my ($self,$x,$y) = @_;

  $self->reset_line();

  my $nx = @$x;
  my $ny = @$y;

  die "*** set_line: Coordinates of different length: $nx $ny\n" if $nx != $ny;

  $self->{n} = $nx;

  my @xnew = @$x;
  my @ynew = @$y;

  $self->{x} = \@xnew;
  $self->{y} = \@ynew;

  $self->{setup} = 1;

  $self->make_unique();		#eliminating double points
  $self->set_xy_min_max();	#sets min/max points
  $self->set_center_point();	#sets center point
  $self->set_area();		#computes area
  my $area = $self->get_area();
  #print STDERR "area: $area\n";
  if( $area < 0 and $self->{must_invert} ) {
    print "inverting line... area negative ($area)\n";
    $self->invert_line();
    $self->set_area();		#computes area
  }
  $self->set_closed();		#sets closed flag, if closed pops last point
  $self->set_convex();		#sets convex flag
  $self->set_latlon();		#sets convex flag
}

sub print_info {

  my ($self,$text,$verbose) = @_;

  print STDERR "$text\n";
  print STDERR "n: $self->{n}   closed: $self->{closed}" .
    "   convex: $self->{convex}\n";
  print STDERR "latlon: $self->{latlon}\n";
  print STDERR "xy-min/max: $self->{xmin} $self->{ymin}" .
    " $self->{xmax} $self->{ymax}\n";
  print STDERR "area: $self->{area}\n";
  print STDERR "center point: $self->{xc} $self->{yc}\n";
}

#--------------------------------------------------------------------

sub is_convex {
  my $self = shift;
  die "*** convex has not been set\n" unless $self->{setup};
  return $self->{convex};
}

sub set_convex {

  my ($self) = @_;

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

#--------------------------------------------------------------------

sub is_closed {
  my $self = shift;
  die "*** closed has not been set\n" unless $self->{setup};
  return $self->{closed};
}

sub set_closed {

  my ($self) = @_;

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

sub make_closed {
  my ($self) = @_;
  $self->{closed} = 1;
}

#--------------------------------------------------------------------

sub invert_line {

  my ($self) = @_;

  my @f = ();
  my $f;

  my $x = $self->{x};
  #print STDERR "x before: @$x\n";
  my @x = reverse(@$x);
  $x = \@x;
  #print STDERR "x after: @$x\n";
  $self->{x} = $x;

  my $y = $self->{y};
  #print STDERR "y before: @$y\n";
  my @y = reverse(@$y);
  $y = \@y;
  #print STDERR "y after: @$y\n";
  $self->{y} = $y;

}

#--------------------------------------------------------------------

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

#--------------------------------------------------------------------

sub get_area {
  my $self = shift;
  die "*** area has not been set\n" unless $self->{setup};
  return $self->{area};
}

sub set_area {
 
  my $self = shift;

  my $n = $self->{n}-1;
  my $x = $self->{x};
  my $y = $self->{y};

  my ($xc,$yc) = $self->get_center_point();
  my $area = 0.;
 
  my ($xm,$ym);
  my ($xn,$yn) = ($x->[$n],$y->[$n]);

  for (my $i=0;$i<=$n;$i++) {
    ($xm,$ym) = ($xn,$yn);
    ($xn,$yn) = ($x->[$i],$y->[$i]);

    $area += areat($xm,$ym,$xn,$yn,$xc,$yc);
  }

  $self->{area} = $area;
  return $area;
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

#--------------------------------------------------------------------

sub is_latlon {
  my $self = shift;
  die "*** latlon has not been set\n" unless $self->{setup};
  return $self->{latlon};
}

sub set_latlon {
  my $self = shift;

  $self->{latlon} = 0;

  if( $self->{xmin} >= -90 and $self->{ymin} >= -90 ) {
    if( $self->{xmax} <= 90 and $self->{ymax} <= 90 ) {
      $self->{latlon} = 1;
    }
  }
}

#--------------------------------------------------------------------

sub get_center_point {
  my $self = shift;
  die "*** center point has not been set\n" unless $self->{setup};
  return ($self->{xc},$self->{yc});
}

sub set_center_point {
 
  my $self = shift;
  my $n = $self->{n};
  my $x = $self->{x};
  my $y = $self->{y};

  my ($xc,$yc) = (0,0);
 
  for( my $i=0;$i<$n;$i++ ) {
    $xc += $$x[$i];
    $yc += $$y[$i];
  }
 
  $self->{xc} = $xc/$n;
  $self->{yc} = $yc/$n;
}

#--------------------------------------------------------------------

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

sub get_xy_min_max {
  my $self = shift;
  die "*** x/y/min/max has not been set\n" unless $self->{setup};
  return ($self->{xmin},$self->{ymin},$self->{xmax},$self->{ymax});
}

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

#--------------------------------------------------------------------

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

  return 0 unless $self->in_min_max($x0,$y0);

  if( $self->{convex} ) {
    return $self->in_convex($x0,$y0);
  } else {
    return $self->in_any_line($x0,$y0);
  }
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

