#!/usr/bin/perl
#
# converts include to modules
#
#-------------------------------------------

$debug = 0;
$do_not_insert = 1;	#for param.h
$do_not_insert = 0;

$what = shift;
$iblock = 0;

init_block();

while(<>) {

  chomp;

  $i++;

  $iimplicit = $i if /^\s*implicit\s*none/;
  $iuse = $i if /^\s*use\s+/;
  $iinclude = $i if /^\s*include\s+\'$what\.h\'/;

  push(@block,$_);

  if( /^\s*end/ ) {
    if( /function/ or /subroutine/ or /module/ or /end\s*$/ ) {
	process_block();
    }
  }
}

process_block();
print STDERR "include blocks found: $iblock\n";

#---------------------------------------------------------

sub process_block {

  print STDERR "processing block: $i\n" if $debug;

  if( $iinclude == -1 ) { 
    write_block();
    init_block();
    return;
  }

  print STDERR "include found: $iinclude\n" if $debug;
  $iblock++;

  if( is_empty($block[$iinclude-1]) and is_empty($block[$iinclude-1]) ) {
    splice(@block,$iinclude,2);
  } else {
    splice(@block,$iinclude,1);
  }

  if( $do_not_insert ) {
    print STDERR "include not substituting...\n";
  } elsif( $iuse != -1 ) {
    splice(@block,$iuse+1,0,"\tuse $what");
  } elsif( $iimplicit ) {
    splice(@block,$iimplicit,0,"\tuse $what\n");
  } else {
    die "cannot parse file\n";
  }

  write_block();
  init_block();
}

#---------------------------------------------------------

sub is_empty {
  my $line = shift;
  if( $line =~ /^\s*$/ ) {
    return 1;
  } else {
    return 0;
  }
}

sub init_block {

  @block = ();

  $i = -1;
  $iimplicit = -1;
  $iuse = -1;
  $iinclude = -1;
}

sub write_block {

  foreach my $line ( @block ) {
    print "$line\n";
  }
}

#---------------------------------------------------------

