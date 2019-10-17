#!/usr/bin/perl
#
# inserts new directory into configuration file
#
# this line must be inserted after if else fi, or shyfem does not compile
# it only inserts once, if run a second time it does not do anything
#
#-----------------------------------------------------

$in_section = 0;

while(<>) {

  if( /^(\s*)\# list files/ ) {
    $in_section = 1;
    print "$1\# ggu list files\n";
  } elsif( $in_section and /^(\s*)fi$/ ) {
    $in_section = 0;
    print;
    print "$1";
    print 'find ${BFMDIR}/src/shyfem -name "*.?90"  -print >> BFM.lst';
    print "\n";
  } else {
    print;
  }
}

#-----------------------------------------------------

