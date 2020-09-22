#!/usr/bin/perl -i.bak
#
# adjusts space after \   
#
#-------------------------------------------

while(<>) {

  #if( /\\ +$/ ) {
  #  print $_;
  #}
  s/\\ +$/\\/;

  print;
}

