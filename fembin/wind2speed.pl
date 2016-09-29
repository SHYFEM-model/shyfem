#!/usr/bin/perl
#
# converts wind file (wx,wy) to speed/dir
#
#----------------------------------------------

$pi = 4*atan2(1.,1.);
$rad = 180 / $pi;

#----------------------------------------------

while(<>) {

  chomp;
  s/^\s+//;
  unless( /^\#/ ) {
    my @f = split;
    my $wx = $f[1];
    my $wy = $f[2];
    my ($s,$d) = convert2polar($wx,$wy);
    $f[1] = $s;
    $f[2] = $d;
    $_ = join(" ",@f);
  }
  print "$_\n";
}

#----------------------------------------------

sub convert2polar
{
  my ($wx,$wy) = @_;

  my ($s,$a,$d);

  $s = sqrt( $wx*$wx + $wy*$wy );

  if( $wx == 0. ) {
    if( $wy > 0. ) {
      $a = 90.;
    } else {
      $a = -90.;
    }
  } else {
    $a = $rad * atan2($wy,$wx);
  }
  $a += 180 if $wx < 0;

  $d = 270 - $a;
  $d -= 360 if $d > 360;
  $d += 360 if $d < 0;

  return ($s,$d);
}

#----------------------------------------------

