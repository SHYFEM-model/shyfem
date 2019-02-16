#!/usr/bin/perl -ws
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# makes structure of program
#
# command line options: 
#
#	-main=routine		reconstruct structure from this routine
#	-maxlevel=#		use this level as maximum descend
#	-upper			return all routines in uppercase
#
#########################################################

%::ignore = ();

%::ignore_ht = (
		 "DEFMAK"		=> 1
		,"MKNAME"		=> 1
		,"DTS2DT"		=> 1
		,"CONWRITE"		=> 1
		,"EXFINTP"		=> 1
		,"EXFFIL"		=> 1
		,"LAGRANGE"		=> 1
		,"CUSTOM"		=> 1
		,"ECOLOGICAL_MODULE"	=> 1
		,"SUBWAVES"		=> 1
		,"SEDI"			=> 1
		,"TURB_CLOSURE"		=> 1
		,"SETDIM"		=> 1

		,"GETPAR"		=> 1
		,"ADDPAR"		=> 1
		,"PUTPAR"		=> 1
		,"SCTPAR"		=> 1
		,"GETFNM"		=> 1
		,"ADDFNM"		=> 1
		,"PUTFNM"		=> 1
		,"SCTFNM"		=> 1

		,"NRDPAR"		=> 1
		,"NLSH2D"		=> 1
		,"NRDLIN"		=> 1

		,"ISCAN"		=> 1
		,"TABULA"		=> 1
		,"ICHAFS"		=> 1
		,"ICHANM"		=> 1
		,"ICHANM0"		=> 1
		,"UPLOW"		=> 1
		,"IROUND"		=> 1

		,"CHKIBC"		=> 1
		,"IGET_BND_ID"		=> 1
		,"BNDSETGET"		=> 1
		,"GET_BND_NAME"		=> 1
		,"GET_BND_IPAR"		=> 1
		,"GET_BND_PAR"		=> 1

		,"IFILEO"		=> 1
		,"IFEMOP"		=> 1
		,"DEF_MAKE"		=> 1
		,"FILNA"		=> 1

		,"IS_R_NAN"		=> 1
		,"CHECK1DR"		=> 1
		,"CHECK2DR"		=> 1

		,"TWB2RH"		=> 1
		,"GASDEV"		=> 1
		,"RGF_SET_GEO"		=> 1
		,"RGF_GET_GEO"		=> 1
		,"SET_SEMI_LAGRANGE"	=> 1
		,"USEUNIT"		=> 1
	  );

%::ignore = %::ignore_ht;

#########################################################

use strict;

$::main = $::main if $::main;
$::maxlevel = 0 unless $::maxlevel;
$::upper = 1 if $::upper;		# per default use lowercase
$::mincount = 0 unless $::mincount;	# min calling count for output

%::defined = ();		# how often is this routine defined
%::functions = ();		# how often is this function defined
%::called = ();			# how often is this routine called
%::calling = ();		# routine calls these routines
%::calledby = ();		# routine is called by these routines
%::routines = ();		# complete info on routines

$::code = "";
$::insection = "";
$::line = "";
$::file = "";
$::space = "";
$::level = 0;
$::maxl = 0;			# maximum level reached
$::maxlevel_reached = 0;

init_ignore();

print "\n===================================================\n";
print "Structure of program:";
print "\n===================================================\n\n";

read_files();
#check_uniquness();
#show_functions();
parse_routines();

print "\n===================================================\n";
print "Statistics of routines found:";
print "\n===================================================\n\n";

foreach my $name ( sort (keys %::defined) ) {

  my $def = $::defined{$name};
  my $called = $::called{$name};
  my $calling = $::calling{$name};
  my $calledby = $::calledby{$name};
  my $item = $::routines{$name};
  my $type = $item->{type};

  $def = defined($def) ? $def : 0;
  $called = defined($called)  ? $called  : 0;
  print "$name   defined: $def   called: $called  type: $type\n";
  my $list1 = make_unique($::calling{$name});
  print "   is calling:   $list1\n";
  my $list2 = make_unique($::calledby{$name});
  print "   is called by: $list2\n";
}

print "\n===================================================\n";
print "Programs called but not defined:";
print "\n===================================================\n\n";

foreach my $name ( sort (keys %::called) ) {
  my $def = exists($::defined{$name}) ? $::defined{$name} : 0;
  my $cal = exists($::called{$name})  ? $::called{$name}  : 0;
  #print "$name     defined: $def      called: $cal\n" if $def == 0;
  print "defined: $def    called: $cal    $name\n" if $def == 0;
}

print "\n===================================================\n";
print "Programs defined but not called:";
print "\n===================================================\n\n";

foreach my $name ( sort (keys %::defined) ) {
  my $def = exists($::defined{$name}) ? $::defined{$name} : 0;
  my $cal = exists($::called{$name})  ? $::called{$name}  : 0;
  #print "$name     defined: $def      called: $cal\n" if $cal == 0;
  print "defined: $def    called: $cal    $name\n" if $cal == 0;
}

#exit 0;

print "\n===================================================\n";
print "Calling sequence:";
print "\n===================================================\n\n";

if( $::main ) {		# only check this routine
  my $name = upperlow($::main);
  sequence($name);
} else {		# do for all not defined routines
  foreach my $name ( sort (keys %::defined) ) {
    my $def = exists($::defined{$name}) ? $::defined{$name} : 0;
    my $cal = exists($::called{$name})  ? $::called{$name}  : 0;
    if( $cal == 0 ) {
      sequence($name);
    }
  }
}

if( $::maxlevel_reached ) {
  print STDERR "max level $::maxlevel has been reached $::maxlevel_reached\n";
}

print_calling_times();

print STDERR "maximum level reached: $::maxl\n";

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub read_files {

 while (<>) {

  if( $::file ne $ARGV ) {
    print STDERR "reading $ARGV\n";
    $::file = $ARGV;
    $::line = 0;
  }

  chomp;

  $::line++;

  my $orig = $_;
  s/^\s*\d+\s+/ /;	# get rid of label number
  s/\s+//g;		# no space
  s/\!.*//g;		# line comment

  if( /^END$/i ) { insert("E",""); }
  if( /^END(FUNCTION|SUBROUTINE)(\w*)/i ) { insert("E",$2); }
  if( /^PROGRAM(\w+)/i ) { insert("P",$1); }
  if( /^SUBROUTINE(\w+)/i ) { insert("S",$1); }
  if( /^FUNCTION(\w+)/i ) { insert("F",$1); }
  if( /^BLOCKDATA(\w+)/i ) { insert("B",$1); }
  if( /^(REAL|INTEGER|DOUBLEPRECISION|LOGICAL)FUNCTION(\w+)/i ) { 
    insert("F",$2); 
  }

  push(@$::code,$orig) if $::code;;
 }
}

sub insert {

  my $type = shift;
  my $name = shift;

  $name = upperlow($name);

  if( $type eq 'E' ) {		# close section
    unless($::insection) {
      print STDERR "file: $::file   line: $::line\n";
      die "end statement found but in no section...\n";
    }
    $name = $::insection;
    my $item = $::routines{$name};
    $item->{code} = $::code;
    $::code = "";
    @::lines = ();
    $::insection = "";
  } else {			# open new section
    if($::insection) {
      print STDERR "file: $::file   line: $::line\n";
      die "routine found but already in section: $::insection\n";
    }
    my @code = ();
    $::code = \@code;
    my $item = {};
    $item->{name} = $name;
    $item->{type} = $type;
    $::routines{$name} = $item;
    $::defined{$name}++;
    $::functions{$name}++ if $type eq "F";
    $::insection = $name;
  }
}

sub handle {

  my $type = shift;
  my $name = shift;

  $name = upperlow($name);

  #print "    ........... call from $::insection to $name ($type)\n";

  $::called{$name}++;
  $::calling{$::insection} .= "$name,";
  $::calledby{$name} .= "$::insection,";
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub show_functions {

  print STDERR "list of functions:\n";
  foreach my $name (keys %::functions) {
    my $count = $::defined{$name};
    print STDERR "$count   $name\n";
  }
}

sub check_uniquness {

  foreach my $name (keys %::defined) {
    my $count = $::defined{$name};
    if( $count != 1 ) {
      print STDERR "multiple definition of routine: $count  $name\n";
    }
  }
}

sub parse_routines {

  foreach my $name (keys %::defined) {
    my $item = $::routines{$name};
    parse_routine($item);
  }
}

sub parse_routine {

  my $item = shift;

  my $name = $item->{name};
  my $type = $item->{type};
  my $code = $item->{code};

  #print STDERR "parsing $name ($type)\n";

  unless( $code ) {
    print STDERR "no code for routine $name\n";
    return;
  }

  $::insection = $name;
  my $nline = 0;

  foreach my $line (@$code) {
    $nline++;
    next if $nline == 1;	# skip definition
    $line =~ s/^\s*\d+\s+/ /;	# get rid of label number

    next if $line =~ /^[cC!]/;		# get rid of line comment

    if( $line =~ /^\s+CALL\s+(\w+)/i ) { handle("C",$1); }
    if( $line =~ /^\s+IF\s*\(.+\)\s*CALL\s+(\w+)/i ) { handle("C",$1); }
    while( $line =~ /(\w+)\s*\(/ig ) { 
      my $func = upperlow($1);
      if( $::functions{$func} ) { 
	#print STDERR "function call: $func in routine $name\n";
        handle("F",$func);
      }
    }

  }
}

#----------------------------------------------------------------

sub sequence {

  my $name = shift;

  return if( $::ignore{$name} );			#programs to ignore

  push_calling_stack($name);
  $::calling_times{$name}++;

  my $def = $::defined{$name};
  my $is_defined = "";
  $is_defined = "(-)" unless $def;

  if( $::maxlevel and $::level > $::maxlevel ) {
    print STDERR "******* maxlevel: $name\n";
    $::maxlevel_reached++;
    #print_calling_stack($name); exit 1;
    return;
  }
  $::level++;

  $::maxl = $::level if $::maxl < $::level;

    print "$::space $name $is_defined\n";
    my $oldspace = $::space;
    $::space .= "   ";
    my $list = $::calling{$name};
    if( $list ) {
      $list =~ s/,$//;
      my @list = split(",",$list);
      foreach my $sub (@list) {
        sequence($sub);
      }
    }
    $::space = $oldspace;

  pop_calling_stack($name);

  $::level--;
}

#----------------------------------------------------------------

sub make_unique {

  my $line = shift;

  return "" unless $line;

  $_ = $line;
  s/,/ /g;
  s/^\s+//g;
  s/\s+/ /g;

  my @f = split;
  my @newlist = ();
  my $equal;

  foreach my $name (@f) {
    $equal = 0;
    foreach my $new (@newlist) {
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

#----------------------------------------------------------------

sub upperlow {
  if( $::upper ) {
    return uc($_[0]);
  } else {
    return lc($_[0]);
  }
}

sub init_ignore {

  foreach (keys %::ignore) {
    my $name = upperlow($_);
    $::ignore{$name} = 1;
  }
}

sub print_calling_times {

  return if $::mincount < 1;

  my @sorted = sort by_count keys %::calling_times;
  #my @sorted = sort keys %::calling_times;

  my $i = 0;
  foreach my $key (@sorted) {
    my $count = $::calling_times{$key};
    next if $count < $::mincount;
    print "calling times (over $::mincount):\n" if $i == 0;
    print "$key   $count\n";
    $i++;
  }
}

sub by_count {

  if( $::calling_times{$a} < $::calling_times{$b} ) {
    return -1;
  } elsif( $::calling_times{$a} > $::calling_times{$b} ) {
    return +1;
  } else {
    return 0;
  }
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub push_calling_stack {

  my $name = shift;

  if( $::calling_stack{$name} ) {
    print_calling_stack();
    die "recursive calling stack: $name\n";
  } else {
    $::calling_stack{$name} = 1;
    push(@::calling_stack,$name);
  }
}

sub pop_calling_stack {

  my $name = shift;

  if( $::calling_stack{$name} ) {
    $::calling_stack{$name} = 0;
    pop(@::calling_stack);
  } else {
    print_calling_stack();
    die "corrupt calling stack: $name\n";
  }
}

sub print_calling_stack {

    my $line = join(" ",@::calling_stack);
    print STDERR "$line\n";
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

