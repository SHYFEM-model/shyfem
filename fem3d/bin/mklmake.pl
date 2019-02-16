#!/usr/bin/perl -wp
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# transforms unix makefile in lahey makefile


  s/\.o\b/.obj/g;

  $exes=1 if /--- EXES\n/;
  $exes=0 if /--- EXES END\n/;

  if( $exes ) {
    if( /^(\w+):/ ) {	#dependency line
	print "$1: $1.\$\(EXE\)\n\n";		#add new rule
	s/^(\w+):/$1.\$\(EXE\):/;		#substitute old
    }
    if( /^\t+\$\(LINKER\)/ ) {	#line for linker
	s/LINKER/LLINKER/;
	s/\$\(LFLAGS\)//;
	s/ -o +\$\@//;
	chop;
	$_ .= ' $(LDEXTRA)' . "\n";
    }
  }


