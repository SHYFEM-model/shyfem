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
# utilities for time series
#
# useage:
#
# require "ts_util.pl";		# ts_util.pl must be in same dir
#
#-------------------------------------

use strict;

$::ts_version = "1.1";
$::flag = -99;		# everything smaller than this is flag

#----------------------------------------------
#
# ts_substitute		substitutes missing values from one to another TS
# ts_interpolate	interpolates missing values in TS
# ts_integrate		integrates TS
# ts_compute_rms	computes rms between two TS
# ts_scale		scales time series like: a = add + mult * a
# 
# ts_get		reads TS
# ts_write		writes TS
# 
# ts_format_value
# ts_set_flag
# ts_is_valid
# ts_version
# 
#----------------------------------------------

sub ts_substitute {

  my ($to,$from) = @_;

  my $n = @$to - 1;

  for my $i (0..$n) {
    my $val = $$to[$i];
    unless( ts_is_valid($val) ) {
      my $newval = $$from[$i];
      #print STDERR "substituting: $i  $val -> $newval\n";
      $$to[$i] = $newval;;
    }
  }
}

sub ts_interpolate {

  my $ts = shift;

  my $n = @$ts;
  my $ivalid = -1;
  my $inblock = 0;

  for(my $i=0;$i<$n;$i++) {
    if( ts_is_valid($$ts[$i]) ) {  #valid
      if( $inblock ) {          #interpolate (end of invalid block found)
        my $istart = $ivalid;
        my $iend = $i;
        ts_interpolate_block($ts,$istart,$iend);
        $inblock = 0;
      }
      $ivalid = $i;
    } else {
     $inblock++;
    }
  }

  if( $inblock ) {              #we are still in an invalid block -> finish
    if( $ivalid == -1 ) {
      die "Cannot interpolate: no good value\n";
    }
    ts_interpolate_block($ts,$ivalid,$n);
  }
}

sub ts_interpolate_block {	# internal subroutine

  my ($ts,$i0,$i1) = @_;

  my $n = @$ts;
  my ($v0,$v1);

  if( $i0 < 0 ) {
    $v0 = $$ts[$i1];
  } else {
    $v0 = $$ts[$i0];
  }
  if( $i1 >= $n ) {
    $v1 = $$ts[$i0];
  } else {
    $v1 = $$ts[$i1];
  }

  my $aux = ($v1-$v0) / ($i1-$i0);

  for( my $i=$i0+1 ; $i<$i1 ; $i++ ) {
    my $val = $$ts[$i];
    my $newval = $v0 + ($i-$i0) * $aux;
    #print STDERR "interpolating $i  $val -> $newval\n";
    $$ts[$i] = $newval;
  }
}

#----------------------------------------------

sub ts_integrate {		

  my ($time,$val,$tstart,$tend) = @_;

  my $n = @$time;
  my $aold = $val->[0];
  my $told = $time->[0];

  $tstart = $time->[0] if not defined $tstart;
  $tend = $time->[$n-1] if not defined $tend;

  my ($apos,$aneg);

  for( my $i=1; $i<$n; $i++ ) {
    my $a = $val->[$i];
    my $t = $time->[$i];
    if( $t > $tstart ) {		# in time window

     if( $told < $tstart ) {
       $aold = ( $a*($tstart-$told) + $aold*($t-$tstart) ) / ($t-$told);
       $told = $tstart;
     }
     if( $t > $tend ) {
       $a = ( $a*($tend-$told) + $aold*($t-$tend) ) / ($t-$told);
       $t = $tend;
     }

     if( $a >= 0 and $aold >= 0 ) {
      my $int = ($t-$told) * 0.5 * ($a+$aold);
      $apos += $int;
     } elsif( $a <= 0 and $aold <= 0 ) {
      my $int = ($t-$told) * 0.5 * ($a+$aold);
      $aneg += $int;
     } else {		# different sign
      my $t0 = ($told*$a+$t*$aold) / ($a+$aold);
      if( $a > 0 ) {
        $aneg -= 0.5 * ($t0-$told) * $aold;
        $apos += 0.5 * ($t-$t0) * $a;
      } else {
        $aneg += 0.5 * ($t0-$told) * $aold;
        $apos -= 0.5 * ($t-$t0) * $a;
      }
     }

    }
    $told = $t;
    $aold = $a;
    last if $told >= $tend;		# finished computing
  }

  my $atot = $apos+$aneg;

  return ($atot,$apos,$aneg);
}

sub ts_compute_rms {

  my ($tres,$zres,$tobs,$zobs,$dt,$tstart,$tend) = @_;

  my $nn = 0;
  my $ddz = 0;
  my $dd2z = 0;

  my $use_interval = $tend - $tstart;
  $use_interval = 0 if $use_interval < 0;	# != 0 only if end > start

  my $nres = @$tres;
  my $nobs = @$tobs;
  my $ires = 0;
  my $iobs = 0;
  my $tr = $tres->[$ires];
  my $to = $tobs->[$iobs];

  while( $ires < $nres and $iobs < $nobs ) {

    if( $tr < $to ) {
      $ires++;
      last if $ires >= $nres;
      $tr = $tres->[$ires];
    } else {
      $iobs++;
      last if $iobs >= $nobs;
      $to = $tobs->[$iobs];
    }

    if( not $use_interval or $tr >= $tstart and $tr <= $tend ) {
      if( abs($to-$tr) < $dt ) {
        $nn++;
        my $zr = $zres->[$ires];
        my $zo = $zobs->[$iobs];
        my $dz = $zr - $zo;
        #print "use: $tr $to   $dz\n";
        $ddz += $dz;
        $dd2z += $dz*$dz;
      }
    }
    last if $use_interval and $tr > $tend;
  }

  if( $nn > 0 ) {
    $ddz /= $nn;
    $dd2z /= $nn;
  }
  my $rms = sqrt( $dd2z );

  return ($nn,$ddz,$rms);
}

#---------------------------------------------------------------

sub ts_scale {

  # scales time series like: a = add + mult * a

  my ($ts,$add,$mult) = @_;

  foreach my $val (@$ts) {
    $val = $add + $mult * $val;
  }

  #my $n = @$ts;
  #for( my $i=0; $i<$n; $i++ ) {
  #}
}

#---------------------------------------------------------------

sub ts_get {		# is not necessary -> use ts_read

  my $file = shift;

  my @time;
  my @val;

  open(FILE,"<$file") || die "*** cannot open file $file\n";

  while(<FILE>) {
    chomp;
    s/^\s+//;
    my ($time,$val) = split;
    push(@time,$time);
    push(@val,$val);
  }

  close(FILE);

  return (\@time,\@val);
}

sub ts_read {

  # reads a variable number of columns

  my $file = shift;

  my $nc = 0;
  my @arrays = ();

  open(FILE,"<$file") || die "*** cannot open file $file\n";

  while(<FILE>) {
    chomp;
    s/^\s+//;
    my @f = split;
    my $m = @f;

    if( $nc == 0 ) {    # allocate arrays
      $nc = $m;
      for( my $i=0;$i<$nc;$i++ ) {
        $arrays[$i] = [];
      }
    } elsif( $nc != $m ) {
      die "*** changing number of colums: $nc $m\n";
    }

    for( my $i=0;$i<$nc;$i++ ) {
      my $a = $arrays[$i];
      push(@$a,$f[$i]);
    }
  }

  close(FILE);

  return @arrays;
}

sub ts_write {

  my @arrays = @_;

  my $nc = @arrays;
  my $ar = $arrays[0];         # first column
  my $n = @$ar;

  for(my $j=0;$j<$n;$j++) {     # loop on length of time series
    my $line = "";
    for( my $i=0;$i<$nc;$i++ ) {
      my $a = $arrays[$i];
      my $val = $a->[$j];
      $line .= "$val ";
    }
    print "$line\n";
  }
}

#----------------------------------------------------------

sub ts_format_value {

  my ($val,$idec,$width) = @_;

  $width = 8 unless $width;

  my $format = "%" . $width . "." . $idec . "f";

  return sprintf($format,$val);
}

sub ts_set_flag {

    $::flag = $_[0] if defined $_[0];

    return $::flag;
}

sub ts_is_valid {

    if( $_[0] > $::flag ) {
      return 1;
    } else {
      return 0;
    }
}

sub ts_version {

  return $::ts_version;
}

#######################################################
1;
#######################################################

