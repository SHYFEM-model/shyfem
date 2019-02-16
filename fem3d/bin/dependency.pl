#!/usr/bin/perl -s
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# -f77   -ifort   -gfortran
# -main				do not flag missing main
#
#------------------------------------------

$file = $ARGV;

while(<>) {

  chomp;

  if( $f77 ) {
	f77();
  } elsif( $gfortran ) {
	f77();
  } elsif( $ifort ) {
	ifort();
  } else {
    print STDERR "unknown compiler or not specified\n";
    exit 1
  }

}

@subs = keys %subs;

foreach $prog (@subs) {

  print "  $prog\n";
}

exit 0;

#------------------------------------------------------

sub ifort {

  if( /^.* undefined reference to \`(\w+)\'$/ ) {
    $prog = $1;
    $prog =~ s/_+$//; 
    return if( $main and $prog eq "MAIN" );
    $subs{$prog}++;
  } elsif( /more undefined references/ ) {
    ;
  } elsif( /: In function \`/ ) {
    ;
  } elsif( /This statement function has not been used/ ) {
    ;
  } else {
    print STDERR "dependency.pl: Cannot process: $_\n";
    exit 1
  }
}

sub f77 {

  if( /^.* undefined reference to \`(\w+)\'$/ ) {
    $prog = $1;
    $prog =~ s/_+$//; 
    return if( $main and $prog eq "MAIN" );
    $subs{$prog}++;
  } elsif( /libgfortranbegin/ ) {
    ;
  } elsif( /more undefined references/ ) {
    ;
  } elsif( /: In function \`/ ) {
    ;
  } elsif( /collect2/ ) {
    ;
  } else {
    print STDERR "dependency.pl: Cannot process: $_\n";
    exit 1
  }
}

