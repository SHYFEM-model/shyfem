#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

while(<>) {

  if( /define LOG_INFO/ ) {
    print; print "    \\\n";
    <>;
    next;
  }

  print;
}

