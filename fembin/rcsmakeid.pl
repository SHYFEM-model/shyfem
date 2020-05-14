#!/usr/bin/perl -i.bak

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

# includes $Id$ in header of file

# first test if header is already there -> must be in the first $nlines lines

$comment = "c";		# comment character
$nlines = 3;		# number of lines that are read for Id

$id = 0;
@lines = ();


$i=0;
while( <> ) {
  $file = $ARGV;
  $i++;
  last if $i > $nlines;
  push(@lines,$_);
  if( /\$Id/ ) {
    $id = 1;
    last;
  }
}

if ( !$id && $i <= $nlines ) {
  die "ERROR $file: too few lines read\n"
}

# if no id we add it

unless( $id ) {
  print "$comment\n";
  print "$comment \$Id\$\n";
}

# print already read lines

while( $_ = shift(@lines) ) {
  print;
}

# print rest of lines

while(<>) {
    print;
}

# reset variables

$id = 0;
@lines = ();

# end of routine

