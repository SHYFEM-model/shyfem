#!/usr/bin/perl -ws
#
# shifts legend of gnuplot
#
# -shift=#		left shift legend by this ammount
# -left			put legend to the left
# -down			put legend at the bottom (entries reversed)
#
# -left and -down can be used contemporarily
#
# the single cuves should be seperated by "% Begin plot" and "% End plot"
#
#-----------------------------------------------------------

use strict;

$::shift = 0 unless $::shift;
$::left = 0 unless $::left;
$::down = 0 unless $::down;

my $ic = 100;
my $dy = 140;
my $dx = 84;
my $xleft = 1500.;
my $ydown = 600.;
my $is_point = 0;	#indicates if we are plotting line or point
my $change_last = 0;	#change last value (for points)

my $last_line;

while(<>) {

  $ic++;

  if( /^% End plot/ and $is_point ) {
    $change_last = 1;
  }
  if( /^% Begin plot/ ) {
    $ic = 0;
    $is_point = 0;
  } else {
    if( $ic == 4 and /setrgbcolor/ ) {	#we probably have to plot a point
      $ic--;				#with points there is one line more
      $is_point = 1;
    }
    if( $::shift ) {
      if( $ic == 4 ) {
        $_ = g_shift($_,$::shift);
      } elsif( $ic == 7 and not $is_point ) {
        $_ = g_subst($_,$::shift);
      } elsif( $change_last ) {
        $last_line = g_subst($last_line,$::shift);
      }
    } else {
      if( $::left ) {
        if( $ic == 4 ) {
          $_ = g_subst($_,$xleft);
        } elsif( $ic == 5 ) {
          s/Rshow/Lshow/;		# left show
        } elsif( $ic == 7 and not $is_point ) {
          $_ = g_subst($_,$xleft-$dx);
        } elsif( $ic == 8 and not $is_point ) {
          $_ = "-" . $_;		# put minus in front
        } elsif( $change_last ) {
          $last_line = g_subst($last_line,$xleft-$dx);
        }
      }
    }
    if( $::down ) {
      if( $ic == 4 ) {
        $_ = g_subst($_,$ydown,1);		# substitute y col
      } elsif( $ic == 7 and not $is_point ) {
        $_ = g_subst($_,$ydown,1);		# substitute y col
	$ydown += $dy;			# one row higher
      } elsif( $change_last ) {
        $last_line = g_subst($last_line,$ydown,1);	# substitute y col
      }
    }
  }

  print $last_line if defined $last_line;	#defer writing
  $last_line = $_;
  $change_last = 0;
}
print $last_line if defined $last_line;

#---------------------------------------------------------

sub g_subst {
  my ($line,$subst,$col) = @_;
  $col = 0 unless defined $col;
  chomp;
  my @f = split(/\s+/,$line);
  $f[$col] = $subst;
  $line = join(" ",@f);
  $line .= "  % substituted\n";
  return $line;
}

sub g_shift {
  my ($line,$shift) = @_;
  chomp;
  my @f = split(/\s+/,$line);
  $f[0] -= $shift;
  $line = join(" ",@f);
  $line .= "  % shifted by $shift\n";
  return $line;
}

#---------------------------------------------------------

