#!/usr/bin/perl -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
#------------------------------------------------

use strict;

my %count = ();
my @noext = ();
my $name;
my $ext;

while( my $file = shift )
{

  next if( -d $file );
  next if( $file =~ m"/tmp/" );
  next if( $file =~ m"/arc/" );
  next if( $file =~ m"/.git/" );

  if( $file =~ m"^.*/(.+)\.([^/]+)" ) {
    $name = $1;
    $ext = $2;
  } elsif( $file =~ m"^.*/\..+" ) {
    # nothing... files like .git etc..
  } else {
    $name = $file;
    $ext = "none";
    push(@noext,$name);
  }

  #print "$file - $name - $ext\n";

  $count{$ext}++;
}

print "stats for extension:\n";
foreach my $key (keys %count) {
  my $c = $count{$key};
  print "$key : $c\n";
}

print "files with no extension:\n";
foreach my $name (@noext) {
  #print "$name ";
}
print "\n";

