#!/usr/bin/perl

# gets header of fortran file

while(<>) {

  if( /^[cC\!]/ || /^\s*$/ ) {			# fortran header
	push(@header,$_);
  } elsif( /^ \*/ || /\/\*/ || /^\s*$/ ) {	# c header
	push(@header,$_);
  } else {
	last;
  }
}

# delete empty lines and "c****" etc.

while( $_ = pop(@header) )
{
  unless( /^.\s*$/ || /^.\**$/ || /^\s*$/ ) {
    push(@header,$_);
    push(@header,"c\n");
    last;
  }
}

print @header;
