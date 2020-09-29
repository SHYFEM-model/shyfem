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
# shifts legend of gnuplot
#
# -shift=#		left shift legend by this ammount
# -left			put legend to the left
# -down			put legend at the bottom (entries reversed)
#
# -left and -down can be used contemporarily
#
# the single curves should be seperated by "% Begin plot" and "% End plot"
#
# now works for versions 4.6 and 5.0
# adjusted for 5.2
#
#-----------------------------------------------------------

use strict;

$::debug = 0 unless $::debug;

$::shift = 0 unless $::shift;	# shift this amount left
$::left = 0 unless $::left;	# shift left with default values
$::down = 0 unless $::down;	# shift down with default values

$::help = 0 unless $::help;	
$::h = 0 unless $::h;	

$::dy = 140;		# vertical distance between entries
$::dx = 84;		# horizontal gap between text and marker
$::xleft = 1300.;	# central x coordinate for legend for left-shifting
$::ydown = 600.;	# vertical start coordinate for down-shifting

if( $::help or $::h ) {
  FullUsage(); exit 0;
} elsif( not $ARGV[0] ) {
  Usage(); exit 1;
}

#---------------------------------------------------------

while(<>) {

  if( /^%%Creator: gnuplot (\S+) patchlevel/ ) {
    $::gversion = $1;
    print STDERR "gnuplot version: $::gversion\n" if $::debug;
  }
  if( /^\/Color (\w+) def/ ) {
    $::color = $1;
    print STDERR "Color = $::color\n" if $::debug;
  }

  if( /^% Begin plot/ ) {
    print;
    skip_intro_lines();
    shift_legend();
    shift_marker();
  } else {
    print;
  }
}

#---------------------------------------------------------

sub shift_marker {

  $_ = <>;

  if( $::gversion >= 5.2 ) {	#we have to find a M command
    print; $_ = <>;
    print; $_ = <>;
  }

  if( /M\s*$/ ) {			# it is a line
    shift_line();
  } else {				# it is a point
    shift_last_point();			# we jump to last point and shift
  }
}

sub shift_last_point {	# we have to shift last point

  my $last_line = $_;

  while(<>) {
    if( /^% End plot/ ) {
      my $end_line = $_;
      if( $::shift ) {
        $last_line = g_shift($last_line,$::shift);
      } else {
        if( $::left ) {
          $last_line = g_subst($last_line,$::xleft-$::dx);
        }
        if( $::down ) {
          $last_line = g_subst($last_line,$::ydown,1);	# substitute y col
          $::ydown += $::dy;			# one row higher
        }
      }
      print $last_line;
      print $end_line;
      return;
    }
    print $last_line;
    $last_line = $_;
  }

}

sub shift_line {	# we have to modify the next two lines

  if( $::shift ) {
    $_ = g_shift($_,$::shift);
  } else {
    if( $::left ) {
      $_ = g_subst($_,$::xleft-$::dx);
    }
    if( $::down ) {
      $_ = g_subst($_,$::ydown,1);		# substitute y col
      $::ydown += $::dy;			# one row higher
    }
  }
  print;

  $_ = <>;
  if( $::left ) {
    chomp;
    $_ = "-" . $_ . " % minus added (ggu)\n";	# put minus in front
  }
  print;
}

sub shift_legend {	# we have to modify the next three lines

  $_ = <>;
  print STDERR "shift: $_" if $::debug;
  if( $::shift ) {
    $_ = g_shift($_,$::shift);
  } else {
    if( $::left ) {
      $_ = g_subst($_,$::xleft);
    }
    if( $::down ) {
      $_ = g_subst($_,$::ydown,1);		# substitute y col
    }
  }
  print;

  $_ = <>;
  if( /Rshow$/ ) {
    s/Rshow/Lshow/ if $::left;
  } else {
    print;
    $_ = <>;
    s/MRshow/MLshow/ if $::left;
    print;
    $_ = <>;
  }
  print;

  $_ = <>;
  print;
}

sub skip_intro_lines {

  my $isrgb = 0;
  my $rgb_max = 1;
  $rgb_max = 2 if $::gversion < 5;
  $rgb_max = 2 if $::color eq "false";	#two setrgbcolor lines for bw

  while(<>) {
    print;
    chomp;
    if( /setrgbcolor$/ ) {
      if( /^L/ ) {	#LC3 setrgbcolor
	$isrgb++;
      } else {		#0.00 0.00 0.00 C LCb setrgbcolor
	last;
      }
    }
    last if $isrgb == $rgb_max;
  }
}

#---------------------------------------------------------

sub g_subst {
  my ($line,$subst,$col) = @_;
  $col = 0 unless defined $col;
  chomp;
  my @f = split(/\s+/,$line);
  my $n = @f;
  $col += $n - 3;	# 3rd last + col
  $f[$col] = $subst;
  $line = join(" ",@f);
  $line .= "  % substituted (ggu)\n";
  return $line;
}

sub g_shift {			# shift left
  my ($line,$shift) = @_;
  chomp;
  my @f = split(/\s+/,$line);
  my $n = @f;
  my $col = $n - 3;	# 3rd last + col
  $f[$col] -= $shift;
  $line = join(" ",@f);
  $line .= "  % shifted by $shift (ggu)\n";
  return $line;
}

#---------------------------------------------------------

sub Usage {
  print STDERR "Usage: gp_shift_legend.pl [-h|-help] [-options] ps-file\n";
}

sub FullUsage {
  Usage();
  print STDERR "  -h|-help      this help screen\n";
  print STDERR "  -left         shift left with default values\n";
  print STDERR "  -down         shift down with default values\n";
  print STDERR "  -shift dx     shift left with dx\n";
  print STDERR "  -debug        write debug messages\n";
}

#---------------------------------------------------------

