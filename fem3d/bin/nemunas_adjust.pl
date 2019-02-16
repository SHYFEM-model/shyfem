#!/usr/bin/perl -wsi.bak
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# adjusts some code for old gfortran compiler
#
#--------------------------------------------------

$undo = 0 unless $undo;

#print STDERR "using undo: $undo\n";

while(<>) {

  if( /(.)NEMUNAS_FIX_(\w+)/ ) {
    #print STDERR "line found ($1-$2): $_";
    $com = $1;
    $what = $2;
    $comment = 1;
    if( $undo and $what eq "NEW" ) { $comment = 0; }
    if( not $undo and $what eq "OLD" ) { $comment = 0; }
    if( $comment ) {
      if( /^(\s*)$com/ ) {
        #nothing ... comment already there
      } else {
        s/^(\s*)/$1$com/;
      }
    } else {
      if( /^(\s*)$com/ ) {
        s/^(\s*)$com/$1/;
      } else {
        #nothing ... already uncommented
      }
    }
    #print STDERR "final: $_";
  }

  print;

}

