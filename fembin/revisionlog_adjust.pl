#!/usr/bin/perl

$FEMDIR=$ENV{FEMDIR}?$ENV{FEMDIR}:$ENV{HOME}/fem;
push(@INC,"$FEMDIR/fembin");

require "revision_getdate.pl";

#--------------------------------------------------
# read in revision log
#--------------------------------------------------

$revfile = shift;

open(FILE,"<$revfile");
@rev = <FILE>;
close(FILE);

#--------------------------------------------------
# insert in hash
#--------------------------------------------------

foreach $line (@rev) {
  $rev{$line} = 1;
}

#--------------------------------------------------
# open VERSION file
#--------------------------------------------------

$version = shift;
#print STDERR "opening file $version...\n";
open(VERS,"<$version");

while($line = <VERS>) {

  $date = get_date($line);
  if( $date ) {
    if( $rev{$line} ) {
      $rev{$line} = 0;
      #print STDERR "**double** $line";
    }
  }
}

#print STDERR "closing file $version...\n";
close(VERS);

@new = ();
foreach $line (@rev) {
  if( $rev{$line} ) {
    push(@new,$line);
  }

  #analyse_line($line);
}

#print @new;

$has_date = 0;
@rev = ();
@old = ();

foreach $line (@new) {
 
  if( is_file_name($line) ) {
    if( $has_date ) {
	push(@rev,@old);
    }
    @old = ();
    $has_date = 0;
  } elsif( is_date($line) ) {
    $has_date = 1;
  }

  push(@old,$line);
}

if( $has_date ) {
  push(@rev,@old);
}

$nref = @rev;
if( $nref <= 0 ) {
  print "No revision log.\n";
} elsif( is_empty_comment($rev[-1]) ) {
  pop(@rev);
}

print @rev;



sub analyse_line {

  my $line = shift;

  if( $ret = is_date($line) ) {
    print STDERR "*** date: $ret\n$line"
  } elsif( $ret = is_empty_comment($line) ) {
    print STDERR "*** comment: $ret\n$line"
  } elsif( $ret = is_file_name($line) ) {
    print STDERR "*** filename: $ret\n$line"
  } else {
    print STDERR "*** other:\n$line"
  }
}

