#!/usr/bin/perl -i.bak
#!/usr/bin/perl
#
# uncomments fortran code

while(<>) {

  if( /^c.*/ ) { next; }	#traditional comment (cC)
  if( /^C.*/ ) { next; }

  s/\![^\'\"]*$/\n/;		#trailing comment

  if( /^\s*\n$/ ) { next; }	#blank line

  &typestat;

  if( $define ) {
	print;
  } elsif( $conti ) {
	print "     +  $statement\n";
  } elsif( $label ) {
	print "$label\t$statement\n";
  } else {
	print "\t$statement\n";
  }

#  print;
}

sub typestat {

  $conti = 0;
  $label = 0;
  $define = 0;

  if( /^\s*\#/ ) {			#compiler directive
	$define = 1;
  } elsif( /^     \S\s*(.*)\n/ ) {	#continuation line
	$conti = 1;
	$statement = $1;
  } elsif( /^\s*\t\s*(.*)\n/ ) {	#regular statement
	$statement = $1;
  } elsif( /^      \s*(.*)\n/ ) {	#regular statement
	$statement = $1;
  } elsif( /^\s*(\d+)\s*(.*)\n/ ) {	#label
	$label = $1;
	$statement = $2;
  } else {				#unknown
	$label = -1;
	$statement = $_;
	chomp($statement);
  }

#  if( $conti )     { print "conti: $statement\n"; }
#  if( $label > 0 ) { print "label: $label - $statement\n"; }
#  if( $label < 0 ) { print "error: $statement\n"; }

  die "error: $statement\n" if( $label < 0 );

}
