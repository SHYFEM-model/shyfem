#!/usr/bin/perl

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use strict;
use date;

my $date = new date;

my $file_obs = $ARGV[0];
my $file_sim = $ARGV[1];

my $year0 = 2014;

#--------------------------------------------------

$date->init_year($year0);              #set reference year

#print STDERR "$file_tot $file_new\n";

my ($tobs,$vobs) = read_file($file_obs);
my ($tsim,$vsim) = read_file($file_sim);

my $nobs = @$tobs;
my $nsim = @$tsim;

#print STDERR "$nobs $nsim\n";

my $rms = rms($tobs,$vobs,$tsim,$vsim);

print STDERR "rms = $rms\n";
print "$rms\n";

#------------------------------------------------

sub rms
{
  my ($tobs,$vobs,$tsim,$vsim) = @_;

  my $ts1 = shift @$tsim;
  my $ts2 = shift @$tsim;
  my $vs1 = shift @$vsim;
  my $vs2 = shift @$vsim;
  my $n = 0;
  my $acum = 0;

  foreach my $to (@$tobs) {
    my $vo = shift(@$vobs);
    next if $to < $ts1;			# no sim for observations
    while( $ts2 < $to ) {
      ($ts1,$vs1) = ($ts2,$vs2);
      $ts2 = shift(@$tsim);
      $vs2 = shift(@$vsim);
      last if not defined $ts2;
    }
    last if not defined $ts2;
    if( $ts1 > $to or $to > $ts2 ) {
      print STDERR "*** error in times: $n  $ts1  $to  $ts2\n";
    }
    $n++;
    my $v = $vs1 + ($vs2-$vs1) * ($to-$ts1)/($ts2-$ts1);
    #print STDERR "$n  $vs1  $v  $vs2\n";
    $acum += ($v-$vo)*($v-$vo);
  }
  $acum /= $n;

  return sqrt($acum);
}

#------------------------------------------------

sub write_file
{
  my ($date,$t,$v) = @_;

  my $file = "new.txt";
  my $n = @$t;

  for( my $i=0; $i<$n; $i++ ) {
    my $tt = $t->[$i];
    my $vv = $v->[$i];
    my ($year,$month,$day,$hour,$min,$sec) = $date->convert_from_it($tt);
    my $line = $date->format_time_date($year,$month,$day,$hour,$min,$sec);
    print "$tt  $vv   $line\n";
  }
  
  print STDERR "$n lines written to stdout\n";
}

sub read_file
{
  my $file = shift;

  my $n = 0;
  my @t = ();
  my @v = ();

  open(FILE,"<$file");

  while(<FILE>) {
    chomp;
    s/^\s+//;
    my @f = split;
    push(@t,$f[0]);
    push(@v,$f[1]);
    $n++;
  }

  close(FILE);

  print STDERR "$n lines read from file $file\n";

  return (\@t,\@v);
}

#------------------------------------------------

