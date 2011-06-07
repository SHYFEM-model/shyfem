#!/usr/bin/perl -si.bak
#
# modifies parameters in files created with gp
#
# modifies file inline
#
# options:
#		-bw   -color
#		-width=#
#
#--------------------------------------------------------

while(<>) {

  chomp;

  if( $width and /^\/gnulinewidth/ ) {
    @f = split;
    $f[1] = $width;
    $_ = join(" ",@f);
  } elsif( $bw and /^\/Color/ ) {
    @f = split;
    $f[1] = "false";
    $_ = join(" ",@f);
  } elsif( $color and /^\/Color/ ) {
    @f = split;
    $f[1] = "true";
    $_ = join(" ",@f);
  }

  print "$_\n";
}
