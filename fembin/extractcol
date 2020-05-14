#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# extracts columns
#
# give file name with -f name, otherwise file name is col
# give x column with -x n, otherwise x col is 1, use 0 to generate index
#
#-------------------------------------------------------------

$newfile = "col";
$xcol = 1;

handle_args();

if( $help ) {
  print STDERR "Usage: extractcol [-help|-h] [-f newfile] [-x xcol] file\n";
  print STDERR "    -h|-help  this help screen\n";
  print STDERR "    -f name   use name for file names\n";
  print STDERR "    -x xcol   use xcol as x-column (0 for automatic index)\n";
  exit 0;
} elsif( !@ARGV ) {
  print STDERR "Usage: extractcol [-help|-h] [-options] file\n";
  exit 1;
}

#print STDERR "newfile=$newfile  xcol=$xcol\n";

$irec = 0;

while(<>) {

  chomp;
  s/^\s+//;		#crop leading white space
  next if /^\#.*/;	#comment line -> ignore
  s/#.*//;		#comment -> delete
  #print STDERR "$_\n";

  $irec++;
  @data = split;

  $nfiles = openfiles() unless( $openfile );

  $x = $irec;
  $x = $data[$xcol-1] if $xcol > 0;
  for($i=1;$i<=$nfiles;$i++) {
    next if $i == $xcol;
    $y = $data[$i-1];
    print $i "$x  $y\n";
  }
}

closefiles();

#--------------------------------------------------------------

sub handle_args {

  while( $arg = $ARGV[0] ) {

    if( $arg =~ /^-f/ ) {
	$newfile = $ARGV[1];
	shift @ARGV;
    } elsif( $arg =~ /^-x/ ) {
	$xcol = $ARGV[1];
	shift @ARGV;
    } elsif( $arg =~ /^-h/ ) {
	$help = 1;
    } elsif( $arg =~ /^-/ ) {
	die "Unknown command line argument: $arg\n";
    } else {
	last;
    }

    shift @ARGV;
  }

}

sub openfiles {

  my $n = @data;
  my $file;
  my $i;

  return $n if $n == 0;

  for($i=1;$i<=$n;$i++) {
    next if $i == $xcol;
    $file = "$newfile.$i";
    open($i,">$file");
  }

  $openfile = 1;

  return $n;
}
  
sub closefiles {

  for($i=1;$i<=$nfiles;$i++) {
    next if $i == $xcol;
    close($i);
  }
}
