#!/usr/bin/perl
#
# converts BFM code to compatible version for compilers ifort and gfortran
#
#---------------------------------------------------------------

$change = 0;

while(<>) {

  chomp;

  if( /^\#INCLUDE/ ) {					#gfortran
    make_change();
    s/^\#INCLUDE/\#include/;
  }

  if( /^\#include\s+\'([\w\.]+)\'/ ) {			#gfortran
    make_change();
    s/^\#include\s+\'([\w\.]+)\'/\#include <$1>/;
  }

  if( /^\#IFDEF/ or /^\#ELSE/ or /^\#ENDIF/ ) {		#gfortran
    make_change();
    s/^\#IFDEF/\#ifdef/;
    s/^\#ELSE/\#else/;
    s/^\#ENDIF/\#endif/;
  }

  if( /integer,static/ ) {				#gfortran
    make_change();
    s/integer,static/integer,save/;
  }

  if( /==\s*\.true\./ ) {				#gfortran
    make_change();
    s/==\s*\.true\./ .eqv. .true./;
  }

  if( /var_ave=0/ ) {					#gfortran
    make_change();
    s/var_ave=0/var_ave=.false./;
  }

  if( /^\#ifdef BFM_NOPOINTERS 1/ ) {			#ifort
    make_change();
    s/^\#ifdef BFM_NOPOINTERS 1/\#ifdef BFM_NOPOINTERS/;
  }

  if( /^\#ifdef DEBUG && BFM_GOTM/ ) {			#ifort
    make_change();
    s/^\#ifdef DEBUG && BFM_GOTM/\#if DEBUG && BFM_GOTM/;
  }

  print "$_\n";

}

if( $change ) {
  exit(1)
} else {
  exit(0)
}

sub make_change		# sets change variable and writes to terminal
{
    $change++;
    print STDERR "$ARGV: $_\n";
}
  
