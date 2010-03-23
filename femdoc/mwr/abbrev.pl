#!/usr/bin/perl

while(<>) {

  chomp;

  if( /^%\s*<MATH_ABBREV>/ ) {
    $inabbrev = 1;
  } elsif( /^%\s*<\/MATH_ABBREV>/ ) {
    $inabbrev = 0;
    &handle_abbrev;	#handle last newcommand
  } elsif( /^\s*\\begin\{document\}/ ) {
    $indocument = 1;
  } elsif( /^\s*\\end\{document\}/ ) {
    $indocument = 0;
  } elsif( $inabbrev ) {
    &handle_abbrev;
  } elsif( $indocument ) {
    &handle_document;
  }
}

print "\\documentclass[12pt]{article}\n";
print "\\usepackage{a4}\n";

foreach $name (@name) {
  $comm = $comm{$name};
  $args = $args{$name};
  print "\\newcommand{\\$name} [$args] {$comm}\n";
}

print "\\begin{document}\n";

foreach $name (@name) {
  $comm = $comm{$name};
  $args = $args{$name};
  $n = $commands{$name};
  $n = "0" unless $n;

  print "\\[\n";
  print "{\\tt $name} \\quad [$args] \\quad $n \\quad \\$name";
  for($i=1;$i<=$args;$i++) {
    print "{\\\#$i}";
  }
  print "\n";
  print "\\]\n";
}



print "\\end{document}\n";

foreach $c (keys %commands) {
#  print "commands: $c $commands{$c}\n";
}

##########################################################

sub handle_document {

  s/%.*$//;	#chop comment

  while( /\\[a-zA-Z]+/ ) {
    $c = $&;
    $c =~ s/\\//;
    $commands{$c}++;
    s/\\[a-zA-Z]+/GGU/;
  }
}

sub handle_abbrev {

  s/%.*$//;	#chop comment

  if( /^\s*\\newcommand/ ) {
    $remember = $_;
    &handle_newcommand($oldcommand);
    $oldcommand = $remember;
  } else {
    $oldcommand .= " $_";
  }
}

sub handle_newcommand {

  $_ = shift;

  if( /^\s*\\newcommand\s*\{\s*\\(\w+)\s*\}\s*\[(\d+)\]\s*\{(.+)\}\s*$/ ) {
    #this is a newcommand with arguments
    $name = $1;
    $args = $2;
    $comm = $3;
  } elsif( /^\s*\\newcommand\s*\{\s*\\(\w+)\s*\}\s*\{(.+)\}\s*$/ ) {
    #this is a newcommand without arguments
    $name = $1;
    $comm = $2;
    $args = "0";
  } elsif( /^\s*$/ ) {
  } else {
    die "bogus newcommand: $_\n";
  }

  if( $name ) {
    push(@name,$name);
    $comm{$name} = $comm;
    $args{$name} = $args;
  }

  print "newcommand: $name [$args] $comm\n" if $name && $debug;
}
