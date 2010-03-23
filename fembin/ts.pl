#!/usr/bin/perl

$what = shift;
$cols = shift;

if( $what eq "-transform" ) {
  $add = shift;
  $mult = shift;
  $func = \&transform;
  print STDERR "debug: $what $cols $add $mult\n";
} elsif( $what eq "-scale" ) {
  $min = shift;
  $max = shift;
  $func = \&scale;
} elsif( $what eq "-scaleaverage" ) {
  $func = \&scaleaverage;
} else {
  die "Unknown transformation: $what\n";
}

@cindex = (1);
@cindex = &make_cols($cols);
print STDERR "columns to alter ($cols) :\n";
&print_array(\@cindex);

@ra = ();
$ncol = 0;

&read_cols;

foreach $cindex (@cindex) {
  $ra = $ra[$cindex];
  $ra[$cindex] = &$func($ra);
}

#&print_array($ra);

&write_cols;

#-------------------------------------------------------

sub transform {

  my $ra = shift;
  my @new = ();

  foreach $val (@$ra) {
    $val = $add + $val * $mult;
    push(@new,$val);
  }

  return \@new;
}

sub scaleaverage {

  my $ra = shift;
  my @new = ();

  my $n = @$ra;
  my $aver = 0;

  foreach $val (@$ra) {
    $aver += $val;
  }
  $aver /= $n;

  #print "average: $aver ($n)\n";

  foreach $val (@$ra) {
    $val = $val / $aver;
    push(@new,$val);
  }

  return \@new;
}

sub scale {

  my $ra = shift;
  my @new = ();

  my $lmin = $lmax = $$ra[0];

  foreach $val (@$ra) {
    $lmin = $val if $val < $lmin;
    $lmax = $val if $val > $lmax;
  }

  #print "min/max: $min $max\n";
  #print "lmin/lmax: $lmin $lmax\n";
  die "min/max: $min $max\n" if $min >= $max;
  die "lmin/lmax: $lmin $lmax\n" if $lmin >= $lmax;

  $dm = ($max-$min) / ($lmax-$lmin);

  foreach $val (@$ra) {
    $val = $min + $dm*($val-$lmin) ;
    push(@new,$val);
  }

  return \@new;
}

#-------------------------------------------------------

sub write_cols {

  my $ra = $ra[0];	#first column
  my $nrow = @$ra;	#number of rows

  #print "ncol = $ncol   nrow = $nrow\n";

  for(my $j=0;$j<$nrow;$j++) {
    for(my $i=0;$i<$ncol;$i++) {
      my $ra = $ra[$i];
      my $val = ${$ra}[$j];
      print "$val ";
    }
    print "\n";
  }

}

sub read_cols {

  while(<>) {

    s/^\s+//;
    my @f = split;
    my $n = @f;

    if( $ncol == 0 ) {
      $ncol = $n;
      for(my $i=0;$i<$n;$i++) {
        my @array = ();
        push(@ra,\@array);
      }
    } elsif( $ncol != $n ) {
      die "Number of columns: $ncol $n\n";
    }

    for(my $i=0;$i<$n;$i++) {
      my $ra = $ra[$i];
      push(@$ra,$f[$i]);
      #print "read_cols: $f[$i]\n";
    }
      
  }
}

sub print_array {

  $ra = shift;

  foreach (@$ra) {
    print STDERR "$_\n";
  }
}

#-------------------------------------------------------

sub make_cols {

  my $cols = shift;

  $cols =~ s/,/ /g;
  $cols =~ s/^\s+//g;
  my @cols = split(/\s+/,$cols);

  my @new = ();

  foreach my $col (@cols) {
    if( $col =~ / (\d+) - (\d+) /x ) {
      die "Cannot yet handle this: $col\n";
    } else {
      push(@new,$col);
    }
  }

  return @new;
}

#-------------------------------------------------------

