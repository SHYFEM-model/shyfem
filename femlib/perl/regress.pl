#!/usr/bin/perl -w
#
# tools for doing regression analysis
#
# how to call:
#
# y = a + b*x
# y = b*x
# y = a + b*exp(c*x);
# y = a + b*x**c;
#
# ($a,$b) = lin_regess($rx,$ry);
# ($a,$b) = lin_regess($rx,$ry,1);	#forcing intercept 0 (through origin)
# $r2 = r2($rx,$ry);			#compute r2
# ($a,$b,$c) = exp_regess($rx,$ry);
# ($a,$b,$c) = poly_regess($rx,$ry);
#
#-----------------------------------------------------------------

use strict;

sub lin_regess {

  my ($rx,$ry,$zero_intercept) = @_;	#with $zero_intercept == 1 make a zero

  my $n = check_arrays("lin_regess",$rx,$ry);

  my $xm = aver($rx);
  my $ym = aver($ry);

  if( $zero_intercept ) {
    $xm = 0; $ym = 0;
  }

  my $xy = 0;
  my $xx = 0;
  for( my $i=0; $i<$n; $i++ ) {
    $xy += ($rx->[$i]-$xm)*($ry->[$i]-$ym);
    $xx += ($rx->[$i]-$xm)**2;
  }

  my $b = $xy / $xx;
  my $a = $ym - $b * $xm;

  return ($a,$b);
}

sub exp_regess {

  my ($rx,$ry) = @_;

  my $n = check_arrays("exp_regess",$rx,$ry);

  my $sk = [];
  my $skold = 0;
  $sk->[0] = $skold;

  for( my $i=1; $i<$n; $i++ ) {
    $skold = $skold + 0.5*($ry->[$i]+$ry->[$i-1])*($rx->[$i]-$rx->[$i-1]);
    $sk->[$i] = $skold;
  }
  
  my $rxx = [];
  my $ryy = [];
  for( my $i=0; $i<$n; $i++ ) {
    $rxx->[$i] = ($rx->[$i]-$rx->[0]);
    $ryy->[$i] = ($ry->[$i]-$ry->[0]);
  }

  my $x2 = 0;
  my $xs = 0;
  my $s2 = 0;
  my $xy = 0;
  my $ys = 0;
  for( my $i=0; $i<$n; $i++ ) {
    $x2 += $rxx->[$i]**2;
    $xs += $rxx->[$i]*$sk->[$i];
    $s2 += $sk->[$i]**2;
    $xy += $rxx->[$i]*$ryy->[$i];
    $ys += $ryy->[$i]*$sk->[$i];
  }

  my ($x,$y) = solve_system($x2,$xs,$xs,$s2,$xy,$ys);

  my $a1 = -$x/$y;
  my $c1 = $y;

  my $c2 = $c1;

  my $phik = 0;
  my $phik2 = 0;
  my $yk = 0;
  my $yphik = 0;
  for( my $i=0; $i<$n; $i++ ) {
    my $phi = exp($c2*$rx->[$i]);
    $phik += $phi;
    $phik2 += $phi**2;
    $yk += $ry->[$i];
    $yphik += $ry->[$i]*$phi;;
  }

  my ($a2,$b2) = solve_system($n,$phik,$phik,$phik2,$yk,$yphik);

  return ($a2,$b2,$c2);
}

sub poly_regess {

  my ($rx,$ry) = @_;

  my $n = check_arrays("poly_regress",$rx,$ry);

  my $rxln = [];
  for( my $i=0; $i<$n; $n++ ) {
    $rxln->[$i] = log( $rx->[$i] );
  }

  return exp_regess($rxln,$ry);
}

#-----------------------------------------------------------------
# auxiliary routines
#-----------------------------------------------------------------

sub solve_system {

# system to solve is:
#
#       a b     x       rx
#     (     ) (   ) = (    )
#       c d     y       ry
#
# (x,y) is returned

  my ($a,$b,$c,$d,$rx,$ry) = @_;

  my $det = $a*$d - $b*$c;

  my $x =  $d * $rx - $b * $ry;
  my $y = -$c * $rx + $a * $ry;

  return ($x/$det,$y/$det);
}

sub r2 {

  my ($rx,$ry) = @_;

  my $n = check_arrays("r2",$rx,$ry);

  my $xm = aver($rx);
  my $ym = aver($ry);

  my $rxy = multiply_arrays($rx,$ry);
  my $rx2 = multiply_arrays($rx,$rx);
  my $ry2 = multiply_arrays($ry,$ry);

  my $xym = aver($rxy);
  my $x2m = aver($rx2);
  my $y2m = aver($ry2);

  my $nom = $xym - $xm*$ym;
  my $denom = ($x2m-$xm**2) * ($y2m-$ym**2);

  return $nom*$nom/$denom;
}

sub aver {

  my $ra = shift;

  my $n = @$ra;
  my $am = 0;
  for( my $i=0; $i<$n; $i++ ) {
    $am += $ra->[$i];
  }
  $am /= $n;

  return $am;
}

sub check_arrays {

  my ($what,$rx,$ry) = @_;

  my $nx = @$rx;
  my $ny = @$ry;
  if( $nx != $ny ) {
    die "$what: different length of arrays: $nx $ny\n";
  }

  return $nx;
}

sub multiply_arrays {

  my ($rx,$ry) = @_;

  my $n = check_arrays("multiply_arrays",$rx,$ry);

  my @xy = ();

  for( my $i=0; $i<$n; $i++ ) {
    my $xy = $rx->[$i] * $ry->[$i];
    push(@xy,$xy);
  }

  return \@xy;
}

sub subtract_arrays {

  my ($rx,$ry) = @_;

  my $n = check_arrays("subtract_arrays",$rx,$ry);

  my @xy = ();

  for( my $i=0; $i<$n; $i++ ) {
    my $xy = $rx->[$i] - $ry->[$i];
    push(@xy,$xy);
  }

  return \@xy;
}

#-----------------------------------------------------------------
# testing routines
#-----------------------------------------------------------------

sub print_array {

  my ($file,$rx,$ry) = @_;

  my $n = check_arrays($file,$rx,$ry);

  open(FILE,">$file");
  for( my $i=0; $i<$n; $i++ ) {
    print FILE "$rx->[$i] $ry->[$i]\n";
  }
  close(FILE);
}

sub make_x {

  my $n = shift;

  my @x = ();

  for( my $i=0; $i<$n; $i++ ) {
    my $r = 2*(rand()-0.5);		# r is in [-1,+1]
    push(@x,$r);
  }

  @x = sort { $a <=> $b } @x;

  return \@x;
}
    
sub perturb {

  my ($ry,$perc) = @_;

  my $ryp = [];

  foreach my $y (@$ry) {
    my $r = $perc*2*(rand()-0.5);		# r is in [-perc,+perc]
    my $yp = $y + $r;
    push(@$ryp,$yp);
  }

  return $ryp;
}

sub make_45_line {

  my ($ry1,$ry2) = @_;

  my @a1 = sort { $a <=> $b } @$ry1;
  my @a2 = sort { $a <=> $b } @$ry2;

  my @n1 = ();
  my @n2 = ();

  my $min = $a1[0];
  $min = $a2[0] if $a2[0] < $min;
  push(@n1,$min);
  push(@n2,$min);

  my $max = $a1[-1];
  $max = $a2[-1] if $a2[-1] > $max;
  push(@n1,$max);
  push(@n2,$max);

  #print STDERR "min/max: $min $max\n";

  return (\@n1,\@n2);
}

#-----------------------------------------------------------------

sub make_lin_reg {

  my ($a,$b,$rx,$perc) = @_;

  my $ry = make_lin_func($a,$b,$rx);
  my $ryp = perturb($ry,$perc);

  my ($a1,$b1) = lin_regess($rx,$ryp);
  my $rr2 = r2($rx,$ryp);
  
  my $ryr = make_lin_func($a1,$b1,$rx);
  my ($r1,$r2) = make_45_line($ryp,$ryr);
  my $rd = subtract_arrays($ryp,$ryr);
  my $aver = aver($rd);

  print_array("lin_orig.txt",$rx,$ry);
  print_array("lin_pert.txt",$rx,$ryp);
  print_array("lin_reg.txt",$rx,$ryr);
  print_array("lin_detrended.txt",$rx,$rd);
  print_array("lin_scatt.txt",$ryp,$ryr);
  print_array("lin_line.txt",$r1,$r2);

  print STDERR "$a $b $a1 $b1 $rr2 $aver\n";
}
  
sub make_lin_func {

  my ($a,$b,$rx) = @_;

  my $ry = [];

  foreach my $x (@$rx) {
    my $y = $a + $b * $x;
    push(@$ry,$y);
  }

  return $ry;
}

#-----------------------------------------------------------------

sub make_exp_reg {

  my ($a,$b,$c,$rx,$perc) = @_;

  my $ry = make_exp_func($a,$b,$c,$rx);
  my $ryp = perturb($ry,$perc);

  my ($a2,$b2,$c2) = exp_regess($rx,$ryp);
  
  my $ryr = make_exp_func($a2,$b2,$c2,$rx);
  my ($r1,$r2) = make_45_line($ryp,$ryr);
  my $rd = subtract_arrays($ryp,$ryr);
  my $aver = aver($rd);

  print_array("exp_orig.txt",$rx,$ry);
  print_array("exp_pert.txt",$rx,$ryp);
  print_array("exp_reg.txt",$rx,$ryr);
  print_array("exp_detrended.txt",$rx,$rd);
  print_array("exp_scatt.txt",$ryp,$ryr);
  print_array("exp_line.txt",$r1,$r2);

  print STDERR "$a $b $c $a2 $b2 $c2 $aver\n";
}
  
sub make_exp_func {		 # y = a + b*exp(c*x);

  my ($a,$b,$c,$rx) = @_;

  my $ry = [];

  foreach my $x (@$rx) {
    my $y = $a + $b * exp($c*$x);
    push(@$ry,$y);
  }

  return $ry;
}

#-----------------------------------------------------------------

sub test_regress {

  my $n = 20;

  my $rx = make_x($n);

  make_lin_reg(2,0.5,$rx,0.1);
  make_exp_reg(0.3,0.6,1.7,$rx,0.1);
}

#-----------------------------------------------------------------
# run test routines
#-----------------------------------------------------------------

test_regress() unless caller;

1;

#-----------------------------------------------------------------
# end of routine
#-----------------------------------------------------------------

