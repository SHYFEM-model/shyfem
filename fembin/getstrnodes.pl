#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# extracts node from STR file: extra, flux, bound

$debug = 1;
$debug = 0;

while(<>) {

  if( /^\s*\&(\w+)/ || /^\s*\$(\w+)/ ) {
    $section = lc($1);
    if( $section eq "end" ) {
      $insection = 0;
      $section = "";
    } else {
      $insection = 1;
    }
    print STDERR "section: $1\n" if $debug;
    next;
  }

  $string = "";
  if( $section eq "extra" ) {
    $string = $section;
    @extra = append($_,@extra);
  } elsif( $section eq "flux" ) {
    $string = $section;
    @flux = append($_,@flux);
  } elsif( $section =~ /^bound(\d+)/ ) {
    $nb = $1;
    $string = "bound ($nb)";
  }

  if( $string ) {
    print STDERR "$string: $_" if $debug;
  }
}

#------------------------------------------------------

$type = 0;
foreach (@extra) {
  print "$_  $type\n";
}

#------------------------------------------------------

push(@flux,0);
unshift(@flux,0);

$type = -1;
foreach (@flux) {
  print "$_  $type\n";
}

#------------------------------------------------------

sub append {

  $_ = shift;
  my @nodes = @_;

  s/,/ /g;
  s/^\s+//;

  my @f = split;

  push(@nodes,@f);

  return @nodes;
}
