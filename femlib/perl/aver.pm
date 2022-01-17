#!/usr/bin/perl
#
# averaging routines
#
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------

sub average_timeseries
{
    my ($ra,$kernel) = @_;

    my $n = @$ra;
    my $n1 = $n - 1;

    my $nk = @$kernel;
    my $ic = int($nk/2);

    my @new = ();

    for(my $i=0;$i<$n;$i++) {
      my $low = $i - $ic;
      $low = 0 if $low < 0;
      my $high = $i + $ic;
      $high = $n1 if $high > $n1;

      my $m = 0;
      my $v = 0;
      for(my $j=$low;$j<=$high;$j++) {
        my $x = $j - $i + $ic;
        my $w = $kernel->[$x];
        $m += $w;
        $v += $w * $ra->[$j];
      }
      $v /= $m if $m;
      push(@new,$v);
    }

    return \@new;
}

sub make_uniform_kernel
{
  my $move = shift;

  my $kernel=[];
  my $fact = 1./(2*$move+1);

  for( my $i=0; $i<=2*$move; $i++ ) {
    $kernel->[$i] = $fact;
  }

  return $kernel;
}

sub make_gaussian_kernel
{
  my $gauss = shift;

  my $kernel=[];
  my $igauss = 3*$gauss;
  my $g2 = $gauss*$gauss;
  my $fact = 1 / sqrt( 2 * 3.14159 * $g2 );

  $kernel->[$igauss] = $fact;
  for( my $i=1; $i<=$igauss; $i++ ) {
    my $x = $i;
    my $w = $fact * exp( -($x*$x)/(2*$g2) );
    $kernel->[$igauss-$i] = $w;
    $kernel->[$igauss+$i] = $w;
  }

  return $kernel;
}

sub print_kernel
{
  my $kernel = shift;

  my $nk = @$kernel;
  my $ic = int($nk/2);
  my $wtot = 0;

  print STDERR "kernel size: $nk $ic\n";
  for( my $i=0; $i<$nk; $i++ ) {
    my $x = $i - $ic;
    my $w = $kernel->[$i];
    $wtot += $w;
    print STDERR "kernel: $i $x $w\n";
  }
  print STDERR "kernel weight: $wtot\n";
}

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------

sub find_min_max
{
  my ($time,$value) = @_;

  my $n = @$value;
  my $grow = 0;

  my @minmax = ();

  for( my $i=1; $i<$n-1; $i++ ) {
    my $lm=$value->[$i-1];
    my $lc=$value->[$i];

    if( $lm < $lc ) {
      if( $grow < 0 ) {
        my $tm = average_time($time,-$grow,$i-1);
	push(@minmax,"$tm $lm min");
      }
      $grow = $i;
    } elsif( $lm > $lc ) {
      if( $grow > 0 ) {
        my $tm = average_time($time,$grow,$i-1);
	push(@minmax,"$tm $lm max");
      }
      $grow = -$i;
    }
  }

  return \@minmax;
}

sub find_min_max_in_min_max
{
  # given minmax already found finds the min/max inside the given windows

  my ($time,$value,$minmax) = @_;

  my $n = @$value;
  my $m = @$minmax;

  my @new_minmax = ();
  my @index = ();

# create index which indicates where the min/max are located

  my $i=0;
  for( my $j=0; $j<$m; $j++ ) {
    my ($t,$l,$w) = split(/\s+/,$minmax->[$j]);
    while( ( $time->[$i] cmp $t ) < 0 ) {
      last if $i >= $n-1;
      $i++;
    }
    $index[$j] = $i;
  }

# find min/max

  for( my $j=0; $j<$m; $j++ ) {
    my ($tc,$lc,$wc) = split(/\s+/,$minmax->[$j]);
    my $ic = $index[$j];
    my $is = $index[$j-1]; $is = 0 if $j <= 0;
    my $ie = $n-1; $ie = $index[$j+1]-1 if $j+1 < $m;
    my $iwhat = 1; $iwhat = -1 if $wc eq "min";
    #print STDERR "$j $is $ic $ie\n";
    die "erroneous index: $j $is $ic $ie\n" if $is>$ic or $ic>$ie;
    my ($ex,$iex,$icount,$iex2) = get_extreme($value,$iwhat,$is,$ie);
    my $t = $time->[$iex];
    $t = average_time($time,$iex,$iex2) if $icount > 1;
    push(@new_minmax,"$t $ex $wc   $icount");
  }

  return \@new_minmax;
}

sub get_extreme
{
  my ($ra,$iwhat,$imin,$imax) = @_;

  my $n = @$ra;

  $iwhat = 1 if $iwhat > 0;
  $iwhat = -1 if $iwhat < 0;
  $imin = 0 unless $imin;
  $imax = $n-1 unless $imax;

  my $ex = $iwhat * $ra->[$imin];
  my $iex = $imin;
  my $iex2 = 0;
  my $icount = 1;

  for( my $i=$imin+1; $i<=$imax; $i++ ) {
    my $v = $iwhat * $ra->[$i];
    if( $v > $ex ) {
      $ex = $v;
      $iex = $i;
      $icount = 1;
    } elsif( $v == $ex ) {
      $iex2 = $i;
      $icount++;
    }
  }

  $ex = $iwhat * $ex;

  return ($ex,$iex,$icount,$iex2);
}

sub average_time 
{
  my ($time,$from,$to) = @_;

  if( $from == $to ) {
    return $time->[$to];
  }

  my $t1 = $::date->unformat_abs($time->[$from]);
  my $t2 = $::date->unformat_abs($time->[$to]);
  my $t = 0.5*($t1+$t2);
  my $tform = $::date->format_abs($t);

  return $tform;
}

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------

sub Test_aver
{
  my $k;

  print STDERR "uniform kernel, move=5\n";
  $k = make_uniform_kernel(5);
  print_kernel($k);

  print STDERR "gaussian kernel, gauss=3\n";
  $k = make_gaussian_kernel(3);
  print_kernel($k);
}

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------

if( $0 =~ /aver.pm$/ ) {
  print "running test routine of aver.pm: $0\n";
  Test_aver();
}

#----------------------------------------------------------
1;
#----------------------------------------------------------
