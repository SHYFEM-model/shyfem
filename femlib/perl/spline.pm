#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#=head1 NAME
#
#    Math::Spline  - Cubic Spline Interpolation of data
#
#=head1 SYNOPSIS
#    
#    require Math::Spline;
#    $spline=new Math::Spline(\@x,\@y)
#    $y_interp=$spline->evaluate($x);
#
#    use Math::Spline qw(spline linsearch binsearch);
#    use Math::Derivative qw(Derivative2);
#    @y2=Derivative2(\@x,\@y);
#    $index=binsearch(\@x,$x);
#    $index=linsearch(\@x,$x,$index);
#    $y_interp=spline(\@x,\@y,\@y2,$index,$x);
#
#=head1 DESCRIPTION
#
#This package provides cubic spline interpolation of numeric data. The
#data is passed as references to two arrays containing the x and y
#ordinates. It may be used as an exporter of the numerical functions
#or, more easily as a class module.
#
#The B<Math::Spline> class constructor B<new> takes references to the
#arrays of x and y ordinates of the data. An interpolation is performed
#using the B<evaluate> method, which, when given an x ordinate returns
#the interpolate y ordinate at that value.
#
#The B<spline> function takes as arguments references to the x and y
#ordinate array, a reference to the 2nd derivatives (calculated using
#B<Derivative2>, the low index of the interval in which to interpolate
#and the x ordinate in that interval. Returned is the interpolated y
#ordinate. Two functions are provided to look up the appropriate index
#in the array of x data. For random calls B<binsearch> can be used -
#give a reference to the x ordinates and the x loopup value it returns
#the low index of the interval in the data in which the value
#lies. Where the lookups are strictly in ascending sequence (e.g. if
#interpolating to produce a higher resolution data set to draw a curve)
#the B<linsearch> function may more efficiently be used. It performs
#like B<binsearch>, but requires a third argument being the previous
#index value, which is incremented if necessary.
#
#=head1 NOTE
#
#requires Math::Derivative module
#
#=head1 EXAMPLE
#
#    require Math::Spline;
#    my @x=(1,3,8,10);
#    my @y=(1,2,3,4);						    
#    $spline=new Math::Spline(\@x,\@y);
#    print $spline->evaluate(5)."\n";
#
#produces the output
#
#2.44    						   
#
#=head1 HISTORY
#
#$Log: Spline.pm,v $
#Revision 1.1  1995/12/26 17:28:17  willijar
#Initial revision
#
#
#=head1 BUGS
#
#Bug reports or constructive comments are welcome.
#
#=head1 AUTHOR
#
#John A.R. Williams <J.A.R.Williams@aston.ac.uk>
#
#=head1 SEE ALSO
#
#"Numerical Recipies: The Art of Scientific Computing"
#W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling.
#Cambridge University Press. ISBN 0 521 30811 9.
#
#=cut

require Exporter;
package spline;
@ISA=qw(Exporter);
@EXPORT_OK=qw(linsearch binsearch spline);
use Carp;
use strict;

sub new {
  my $type=shift;
  my $self=[];
  push @{$self},shift; # x
  push @{$self},shift; # y
  my $y2=[Derivative2($self->[0],$self->[1])];
  push @{$self},$y2;
  bless $self,$type;
}

sub evaluate {
  my ($self,$v)=@_;
  my $idx=binsearch($self->[0],$v);
  spline($self->[0],$self->[1],$self->[2],$idx,$v);
}

sub spline { 
  my ($x,$y,$y2,$i,$v)=@_;
  my ($klo,$khi)=($i,$i+1);
  my $h=$x->[$khi]-$x->[$klo];
  if ($h==0) { croak "Zero interval in spline data.\n"; }
  my $a=($x->[$khi]-$v)/$h;
  my $b=($v-$x->[$klo])/$h;
  return $a*$y->[$klo] + $b*$y->[$khi]
      +(($a*$a*$a-$a)*$y2->[$klo]
	+($b*$b*$b-$b)*$y2->[$khi])*($h*$h)/6.0;
}

sub binsearch { # binary search routine finds index just below value
  my ($x,$v)=@_;
  my ($klo,$khi)=(0,$#{$x});
  my $k;
  while (($khi-$klo)>1) {
    $k=int(($khi+$klo)/2);
    if ($x->[$k]>$v) { $khi=$k; } else { $klo=$k; }
  }
  return $klo;
}

sub linsearch { # more efficient if repetatively doint it
  my ($x,$v,$khi)=@_; $khi+=1;
  my $n=$#{$x};
  while($v>$x->[$khi] and $khi<$n) { $khi++; }
  $_[2]=$khi-1;
}

sub Derivative1 {
    my ($x,$y)=@_;
    my @y2;
    my $n=$#{$x};
    $y2[0]=($y->[1]-$y->[0])/($x->[1]-$x->[0]);
    $y2[$n]=($y->[$n]-$y->[$n-1])/($x->[$n]-$x->[$n-1]);
    my $i;
    for($i=1; $i<$n; $i++) {
	$y2[$i]=($y->[$i+1]-$y->[$i-1])/($x->[$i+1]-$x->[$i-1]);
    }
    return @y2;
}

sub Derivative2 {
    my ($x,$y,$yp1,$ypn)=@_;
    my $n=$#{$x};
    my ($i,@y2,@u);
    if (!defined $yp1) {
	$y2[0]=0; $u[0]=0;
    }
    else {
	$y2[0]=-0.5;
	$u[0]=(3/($x->[1]-$x->[0]))*(($y->[1]-$y->[0])/($x->[1]-$x->[0])-$yp1);
    }
    for($i=1; $i<$n; $i++) {
	my $sig=($x->[$i]-$x->[$i-1])/($x->[$i+1]-$x->[$i-1]);
	my $p=$sig*$y2[$i-1]+2.0; 	
	$y2[$i]=($sig-1.0)/$p;
	$u[$i]=(6.0*( ($y->[$i+1]-$y->[$i])/($x->[$i+1]-$x->[$i])-
		      ($y->[$i]-$y->[$i-1])/($x->[$i]-$x->[$i-1])
		     )/
		($x->[$i+1]-$x->[$i-1])-$sig*$u[$i-1])/$p;
    }
    my ($qn,$un);
    if (!defined $ypn) {
	$qn=0;
	$un=0;
    }
    else {
	$qn=0.5;
	$un=(3.0/($x->[$n]-$x->[$n-1]))*
	    ($ypn-($y->[$n]-$y->[$n-1])/($x->[$n]-$x->[$n-1]));
    }
    $y2[$n]=($un-$qn*$u[$n-1])/($qn*$y2[$n-1]+1.0);
    for($i=$n-1; $i>=0; --$i) {
	$y2[$i]=$y2[$i]*$y2[$i+1]+$u[$i];
    }
    return @y2;
}

sub intp_neville {	# has still to be tested

  my ($nintp,$xa,$ya,$x) = @_;

  my $n = $nintp - 1;

  my @p = @$ya;

  for( my $k = 1; $k<=$n; $k++ ) {
    for( my $i = $n; $i>=$k; $i-- ) {
      my $xl = $xa->[$i-$k];
      my $xh = $xa->[$i];
      $p[$i] = ( ($x-$xl)*$p[$i] - ($x-$xh)*$p[$i-1] ) / ($xh-$xl)
    }
  }

  return $p[$n];
}

#------------------------------------------------
1;
#------------------------------------------------

