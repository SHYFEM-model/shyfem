#!/usr/bin/perl -s
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# formats update list of cvs
#
# recognizes options:  -detail  -with_date

while(<>) {

  ($what,$file) = split;

  if( $what eq "M" ) {
	push(@modified,$file);
  } elsif( $what eq "R" ) {
	push(@removed,$file);
  } elsif( $what eq "A" ) {
	push(@added,$file);
  } elsif( $what eq "?" ) {
	push(@unknown,$file);
  } else {
	push(@various,$file);
  }

}

@next = @modified;
if( @next ) {
  print "   Modified files:\n";
  if( $detail ) {
    &format_detailed_list(@next);
  } else {
    &format_list(@next);
  }
}

@next = @removed;
if( @next ) {
  print "   Removed files:\n";
  &format_list(@next);
}

@next = @added;
if( @next ) {
  print "   Added files:\n";
  &format_list(@next);
}

@next = @various;
if( @next ) {
  print "   Unknown mode files:\n";
  &format_list(@next);
}

sub format_detailed_list {

  use FileHandle;

  my $file;
  my $changes;

  my $date = "11.22.3333";

  if( $with_date ) {
    format_name     STDOUT "WITH_DATE";
  }

  while( $file = shift ) {
    $changes = `cvs diff -n $file | tail -n +6 | d+.pl`;
    $date = format_date($file);
    write;
    #print "      $file: $changes";
  }

format STDOUT =
@>>>>>>>>>>>> @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$file,         $changes
.

format NO_DATE =
@>>>>>>>>>>>> @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$file,         $changes
.

format WITH_DATE =
@>>>>>>>>>>>> @>>>>>>>>>> @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$file,        $date,      $changes
.
}

sub format_list {

  my $file;
  my $length = 0;
  my $l;

  print "      ";
  while( $file = shift ) {
    $l = length($file) + 1;
    if( $length + $l > 60 ) { 
	print "\n"; 
	print "      ";
	$length = 0;
    }
    $length += $l;
    print "$file ";
  }
  print "\n";
}

sub format_date {

    my $file = shift;

    my @f = stat($file);
    my @t = localtime($f[9]);

    my $day = "$t[3]";
    $day = "0" . $day if( $day < 10 );
    my $month = "$t[4]";
    $month++;	#month is given from 0-11
    $month = "0" . $month if( $month < 10 );
    my $year = 1900 + $t[5];

    my $time = "$day.$month.$year";
    #print "date: $time\n";

    return $time;
}

