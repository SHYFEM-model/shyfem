#!/usr/bin/perl -ws
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# computes average
#
# -move=m	computes moving average over m data, else arithmetic average
# -col=c	averages only column c, else all columns
# -noxcol	files has no x/time column
# -fact=f	multiplies columns with f
# -sum		computes total and not average
# -std		computes standard deviation
# -min		computes minimum
# -max		computes maximum
# -minmax	computes minimum and maximum
# -averstd	computes average and standard deviation
# -format=val	number of significant digits (e.g., 0.01: 45.34, 1: no fract)
#
#------------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use warnings;
use strict;
use date;

# handle options

$::debug = 0;

$::h = 0 unless $::h;
$::help = 0 unless $::help;
$::move = 0 unless $::move;
$::regress = 0 unless $::regress;
$::noxcol = 0 unless $::noxcol;
$::col = 0 unless $::col;
$::fact = 0 unless $::fact;
$::sum = 0 unless $::sum;
$::std = 0 unless $::std;
$::min = 0 unless $::min;
$::max = 0 unless $::max;
$::minmax = 0 unless $::minmax;
$::averstd = 0 unless $::averstd;
$::format = 0 unless $::format;

my @files = @ARGV;
my $nfiles = @files;

fullusage() if $::h or $::help;
usage() unless $nfiles;

print STDERR "total number of files: $nfiles\n" if $nfiles > 1;

$::date = new date;

my $total = 0;

if( $::move ) {
    my @new = <>;
    my @cols = read_cols(\@new);
    print STDERR "computing moving average with $::move\n";
    $cols[1] = maver($::move,$cols[1]);
    print_cols(@cols);
} elsif( $::regress ) {
    print STDERR "computing linear regression\n";
    my @new = <>;
    my @cols = read_cols(\@new);
    my $n = @cols;
    my $time = convert_date($cols[0]);
    for(my $i=1;$i<$n;$i++) {
      my ($b0,$b1,$df,$t) = regress($time,$cols[$i]);
      print STDERR "col $i:  $b0 $b1 $df $t\n";
      write_regress($i,$time,$b0,$b1);
    }
    print "$::regline0\n";
    print "$::regline1\n";
} else {
    my @new = <>;
    my @cols = read_cols(\@new);
    my $n = @cols;
    $n-- if is_date($cols[$n-1]);
    my ($colmin,$colmax);
    if( $::noxcol ) {		#no time column
      $colmin = 0;
    } else {
      $colmin = 1;
    }
    if( $::col ) {		#col 0 is time column
      $colmin = $::col;
      $colmax = $::col+1;
    } else {
      $colmax = $n;
    }
    for(my $i=$colmin;$i<$colmax;$i++) {
        my ($aver,$std) = aver($cols[$i]);
	if( $::std ) {
          print_format($std);
	} elsif( $::averstd ) {
          print_format($aver); print " +- "; print_format($std);
	} elsif( $::min ) {
          print_format($aver);
	} elsif( $::max ) {
          print_format($std);
	} elsif( $::minmax ) {
          print_format($aver); print " - "; print_format($std);
        } else {
          print_format($aver);
        }
	print "  ";
    }
    print "\n";
}

###############################################################

sub aver
{
    my $ra = shift;

    my $total = 0;
    my $accum = 0;
    my $min = 0;
    my $max = 0;
    my $n = @$ra;

    return (0,0) unless $n;

    $min = $ra->[0];
    $max = $ra->[0];
    foreach (@$ra) {
	$total += $_;
	$min = $_ if $_ < $min;
	$max = $_ if $_ > $max;
    }

    return ($min,$max) if $::min or $::max or $::minmax;

    my $aver = $total / $n;

    foreach (@$ra) {
	$accum += ($_-$aver)*($_-$aver);
    }
    
    my $std = sqrt( $accum / $n );

    $aver = $total if $::sum;

    return ($aver,$std);
}

###############################################################

sub maver
{
    my ($move,$ra) = @_;

    my $n = @$ra;
    my $n1 = $n - 1;

    my @new = ();

    for(my $i=0;$i<$n;$i++) {
      #print STDERR "moving average $i\n" if $i%100 == 0;
      my $low = $i - $move;
      $low = 0 if $low < 0;
      my $high = $i + $move;
      $high = $n1 if $high > $n1;
      my $m = 0;
      my $v = 0;
      #print STDERR "moving average $i ($low,$high)\n" if $i%100 == 0;
      for(my $j=$low;$j<=$high;$j++) {
	$m++;
	$v += $$ra[$j];
      }
      $v /= $m if $m;
      #print STDERR "moving average $i ($low,$high) $m -> $v\n" if $i%100 == 0;
      push(@new,$v);
    }

    return \@new;
}

###############################################################

sub regress
{
  my ($time,$vals) = @_;

  my $n = @$time;
  return (0,0,0,0) if $n == 0;

  my ($xm,$ym) = (0,0);
  my ($sxy,$sxx) = (0,0);

  for(my $i=0; $i<$n; $i++) {
    $xm += $time->[$i];
    $ym += $vals->[$i];
  }
  $xm /= $n;
  $ym /= $n;

  for(my $i=0; $i<$n; $i++) {
    my $dx = $time->[$i] - $xm;
    my $dy = $vals->[$i] - $ym;
    $sxy += $dx*$dy;
    $sxx += $dx*$dx;
  }

  my $b1 = $sxy/$sxx;
  my $b0 = $ym - $b1 * $xm;

  my $ssr = 0;
  for(my $i=0; $i<$n; $i++) {
    my $x = $time->[$i];
    my $y = $vals->[$i];
    my $ye = $b0 + $b1 * $x;
    my $dr = $y - $ye;
    $ssr += $dr * $dr;
  }

  my $df = $n - 2;	#degree of freedom
  my $t = 0;
  if( $df > 0 and $ssr > 0 ) {
    $t = ( $b1 * sqrt($df) ) / sqrt( $ssr / $sxx );
  }

  return($b0,$b1,$df,$t);
}

sub write_regress
{
  my ($i,$time,$b0,$b1) = @_;

  my $t0 = $time->[0];
  my $t1 = $time->[-1];

  my $y0 = $b0 + $b1 * $t0;
  my $y1 = $b0 + $b1 * $t1;

  $y0 = put_format($y0);
  $y1 = put_format($y1);

  $t0 += $::atime0;
  $t1 += $::atime0;

  my $ta0 = $::date->format_abs($t0);
  my $ta1 = $::date->format_abs($t1);
  
  #print "$ta0    $y0\n";
  #print "$ta1    $y1\n";

  if( $i == 1 ) {
    $::regline0 = "$ta0  ";
    $::regline1 = "$ta1  ";
  }
  $::regline0 .= "$y0 ";
  $::regline1 .= "$y1 ";
}

###############################################################

sub print_cols
{
    my @cols = @_;

    my $col = $cols[0];
    my $nrows = @$col;

    for(my $i=0;$i<$nrows;$i++) {
      foreach $col (@cols) {
        my $value = $$col[$i];
        print "$value ";
      }
      print "\n";
    }
}

###############################################################

sub convert_date
{
  my $timecol = shift;

  my $n = @$timecol;

  $::atime0 = -1;
  $_ = $timecol->[0];
  return $timecol unless /[:-]+/;
  $::atime0 = $::date->unformat_abs($_);

  my @time = ();

  for(my $i=0; $i<$n; $i++) {
    $_ = $timecol->[$i];
    my $atime = $::date->unformat_abs($_);
    push(@time,$atime-$::atime0);
  }

  return \@time;
}

###############################################################

sub read_cols
{
    my $lines = shift;

    my @cols = ();

    foreach (@$lines) {
	s/^\s+//;
	next if /^\#/;
	my @f = split;
        my $ncols = @f;
        for(my $i=0;$i<$ncols;$i++) {
          my $ra = $cols[$i];
	  unless( $ra ) {
            my @new = ();
	    $cols[$i] = \@new;
            $ra = $cols[$i];
	  }
	  $f[$i] *= $::fact if $i > 0 and $::fact != 0;
          push(@$ra,$f[$i]);
        }
    }

    return @cols;
}

sub put_format {

  my $val = shift;

  if( $::format ) {
    my $f = $::format;
    $val /= $f;
    if( $val > 0 ) {
      $val = int($val+0.5);
    } else {
      $val = int($val-0.5);
    }
    $val *= $f;
  }

  return "$val";
}

sub print_format {

  my $val = shift;

  $val = put_format($val);

  print "$val";
}

###############################################################

sub is_date
{
  my $ra = shift;

  if( $ra->[0] =~ /::/ ) {
    print STDERR "last column is date column\n" if $::debug;
    return 1;
  } else {
    return 0;
  }
}

###############################################################

sub usage
{
  print "Usage: average [-h|-help] [options] file(s)\n";
  exit 1;
}

sub fullusage
{
  print "Usage: average [-h|-help] [options] file(s)\n";
  print "  options:\n";
  print "  -move=m	computes moving average over m data (on both sides)\n";
  print "  -regress	computes linear regression\n";
  print "  -col=c	averages only column c, else all columns\n";
  print "  -noxcol	file has no x/time column\n";
  print "  -fact=f	multiplies columns with f\n";
  print "  -sum		computes total and not average\n";
  print "  -std		computes standard deviation\n";
  print "  -min		computes minimum\n";
  print "  -max		computes maximum\n";
  print "  -minmax	computes minimum and maximum\n";
  print "  -averstd	computes average and standard deviation\n";
  print "  -format=val	number of significant digits \n";
  print "                  (e.g., 0.01: 45.34, 1: no fract)\n";
  exit 1;
}

###############################################################

