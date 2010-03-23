#!/usr/bin/perl
#
# makes structure of program

#########################################################
%ignore = (
		 "WEXIT"		=> 1
		,"WERR"			=> 1
		,"WMESS"		=> 1
		,"WSCREEN"		=> 1
		,"WLABELS"		=> 1
		,"WUPDATE"		=> 1
		,"GETMES"		=> 1
		,"BRKERR"		=> 1
		,"PROMPT"		=> 1
	  );
%force =  (
		 "WASPB"		=> 1
	  );
#########################################################

print "\nStructure of program:\n\n";

while (<>) {

  if( /^\s+PROGRAM\s+(\w+)/i ) {
    &handle("P",$1);
  }

  if( /^\s+SUBROUTINE\s+(\w+)/i ) {
    &handle("S",$1);
  }

  if( /^\s+FUNCTION\s+(\w+)/i ) {
    &handle("F",$1);
  }

  s/^\s*\d+\s+/ /;	#get rid of label number

  if( /^\s+CALL\s+(\w+)/i ) {
    &handle("C",$1);
  }

  if( /^\s+IF\s*\(.+\)\s*CALL\s+(\w+)/i ) {
    &handle("C",$1);
  }

}

print "\nStatistics of programs found:\n\n";

foreach $name ( sort (keys %mentioned) ) {
  $def = exists($defined{$name}) ? $defined{$name} : 0;
  $cal = exists($called{$name})  ? $called{$name}  : 0;
  print "$name     $def      $cal\n";
  $list = &unique($calling{$name});
  print "   is calling:   $list\n";
  $list = &unique($calledby{$name});
  print "   is called by: $list\n";
}

print "\nPrograms not defined:\n\n";

foreach $name ( sort (keys %mentioned) ) {
  $def = exists($defined{$name}) ? $defined{$name} : 0;
  $cal = exists($called{$name})  ? $called{$name}  : 0;
  print "$name     $def      $cal\n" if $def == 0;
}

print "\nPrograms not called:\n\n";

foreach $name ( sort (keys %mentioned) ) {
  $def = exists($defined{$name}) ? $defined{$name} : 0;
  $cal = exists($called{$name})  ? $called{$name}  : 0;
  print "$name     $def      $cal\n" if $cal == 0;
}

print "\nCalling sequence:\n\n";

foreach $name ( sort (keys %mentioned) ) {
  $def = exists($defined{$name}) ? $defined{$name} : 0;
  $cal = exists($called{$name})  ? $called{$name}  : 0;
  if( $cal == 0 ) {
    &sequence($name);
  }
}

##################################

sub handle {

  my $type = shift;
  my $name = shift;

  $name = uc($name);

  #print "........... $type $name\n";

  if( $type eq "C" ) {
    print "     calls   $name\n";
    $mentioned{$name}++;
    $called{$name}++;
    $calling{$insection} .= "$name,";
    $calledby{$name} .= "$insection,";
  } else {
    print "$type   $name\n";
    $mentioned{$name}++;
    $defined{$name}++;
    $insection = $name;
  }
}

sub sequence {

  my $name = shift;
  my $oldspace;

  if( $ignore{$name} ) {			#programs to ignore
  } elsif( $force{$name} || $defined{$name} ) {	#not progs that are not defined
    print "$space$name\n";
    $oldspace = $space;
    $space .= "   ";
    my $list = $calling{$name};
    $list =~ s/,$//;
    my @list = split(",",$list);
    foreach $sub (@list) {
      &sequence($sub);
    }
    $space = $oldspace;
  }
}

sub unique {

  my $line = shift;

  $_ = $line;
  s/,/ /g;
  s/^\s+//g;
  s/\s+/ /g;

  my @f = split;
  my @newlist = ();
  my $equal;

  foreach $name (@f) {
    $equal = 0;
    foreach $new (@newlist) {
      if( $new eq $name ) {
        $equal = 1;
        last;
      }
    }
    unless( $equal ) {
      push(@newlist,$name)
    }
  }

  $line = join(" ",@newlist);

  return $line;
}
