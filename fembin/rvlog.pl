#!/usr/bin/perl

# gets revision log from header file (ASCII)

while(<>) {

  next if( /^\s*$/ );	#empty line

  if( /^c/ ) {		#comment
    if( /^c revision log :/ ) {
      &new_rev;
    }
  } else {		#file name
    chop;
    $file = $_;
  }
}

sub new_rev {

  @rev = ();

  # append file name on revision log line

  chop;
  s/\s*$//;
  $_ .= " $file\n";

  # save first two lines

  push(@rev,$_);
  $_ = <>;
  push(@rev,$_);
  $_ = <>;

  # loop until break line found

  do {
    push(@rev,$_);
    $_ = <>;
  } until( /^c?\s*\n/ );

  print "\n";
  print @rev;
  print "\n";
}

