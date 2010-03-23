#!/usr/bin/perl
#
# gets name of subroutines/functions defined in file


while(<>) {

	chomp;

	if( /^\s+subroutine\s+(\w+)/i ) {
	  print "$1\n";
	} elsif( /^\s+function\s+(\w+)/i ) {
	  print "$1\n";
	} elsif( /^\s+([\w\s]+)function\s+(\w+)/i ) {
	  $name = $2;
	  $mod = $1;
	  $mod =~ s/\s+//g;
	  if( $mod =~ /^(integer|real|logical|character|doubleprecision)$/i ) {
		print "$name\n";
	  }
	}
}

