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
#-----------------------------------

use strict;

my $home = $ENV{"HOME"};

$::copyright = "$home/shyfem/femcheck/copyright/copyright_notice.txt";
$::copyshort = "$home/shyfem/femcheck/copyright/copyright_short.txt";

$::type = "" unless $::type;

my $file = $ARGV[0];
$file = "" unless $file;

my $type = $::type;
$type = find_file_type($file) unless $type;

if( not $file ) {
  print STDERR "no file given...\n";
  die "Usage: include_copyright.pl file\n";
}
if( not $type ) {
  print STDERR "cannot determine file type of file $file\n";
  die "Usage: include_copyright.pl file\n";
}

my $copy = read_copyright($type);
$copy = adjust_copyright($copy,$type);

print_file($type,$copy);

#-----------------------------------

sub adjust_copyright
{
  my ($copy,$type) = @_;

  my ($char,$bline,$eline,$line);

  if( $type eq "fortran" ) {
    $char = "!";
    my $stars = "-" x 74;
    $bline = "!" . $stars . "\n";
    $eline = $bline;
  } elsif( $type eq "c" ) {
    $char = " *";
    my $stars = "*" x 72;
    $bline = "/" . $stars . "\\\n";
    $eline = "\\" . $stars . "/\n";
  } elsif( $type eq "tex" ) {
    $char = "%";
    my $stars = "-" x 72;
    $bline = "%" . $stars . "\n";
    $eline = $bline;
  } elsif( $type eq "text" ) {
    $char = "#";
    my $stars = "-" x 72;
    $bline = "#" . $stars . "\n";
    $eline = $bline;
  } elsif( $type eq "script" ) {
    $char = "#";
    my $stars = "-" x 72;
    $bline = "#" . $stars . "\n";
    $eline = $bline;
  } else {
    die "unknown type: $type\n";
  }

  my @new = ();

  push(@new,"\n");
  push(@new,$bline);

  foreach (@$copy) {
    $line = $char . $_;
    push(@new,$line);
  }

  push(@new,$eline);
  push(@new,"\n");

  return \@new;
}

sub read_copyright
{
  my $type = shift;

  my $file;

  if( $type eq "fortran" or $type eq "c" or $type eq "tex" ) {
    $file = $::copyright;
  } else {
    $file = $::copyshort;
  }

  open(COPY,"<$file") || die "cannot read file: $file\n";
  my @copy = <COPY>;
  close(COPY);

  return \@copy;
}

sub print_copyright
{
  my $a = shift;

  foreach (@$a) {
    print "$_";
  }
}

sub print_file
{
  my ($type,$copy) = @_;

  my $debug = 0;	# set to > 0 for debug
  my $header = 1;	# we are still in header
  my $oldcopy = 0;	# we are in old copyright

  if( $type eq "script" ) {
    $_ = <>; print;
    $_ = <>;
    handle_2line_script($copy);
  }

  print_copyright($copy);

  while(<>) {
    if( $header ) {
      if( $type eq "fortran" ) {
        next if /^[cC!]\s*$/;
        next if /^\s*$/;
        next if /^[cC!]\s+\$Id:/;
      } elsif( $type eq "c" ) {
	next if /^\/\* \$Id:/;
        $oldcopy = 1 if /Copyright \(c\)/;
        if( /Revision / ) {
          $oldcopy = 0;
	  $header = 0;
        }
	next if $oldcopy;
	print;
	next;
      } elsif( $type eq "tex" ) {
        next if /^\s*$/;
        next if /^%\s+\$Id:/;
      } elsif( $type eq "text" ) {
        next if /^#\s*$/;
        next if /^#\s+\$Id:/;
      } elsif( $type eq "script" ) {
        next if /^#\s*$/;
        next if /^#\s+\$Id:/;
      } else {
	die "unknown type: $type\n";
      }
    }
    last if( --$debug == 0 );	# only executed if debug was > 0
    print;
    $header = 0;
  }
}

sub handle_2line_script
{
  my $copy = shift;

  if( /^#\s*$/ ) {		# empty line with leading #
    $copy->[0] = "#\n";
    $copy->[-1] = "#\n";
  } elsif( not /^\s*$/ ) {	# something on line
    $copy->[-1] .= "$_";
  } else {
    # we simply drop the empty line
  }
}

sub find_file_type
{
  my $file = shift;

  return "" unless $file;

  my $type = "";

  $_ = $file;

  if( /\.f$/ or /\.f90$/ or /\.F90$/ ) {
    $type = "fortran";
  } elsif( /\.c$/ ) {
    $type = "c";
  } elsif( /\.tex$/ ) {
    $type = "tex";
  } elsif( /\.sh$/ or /\.pl$/ or /\.pm$/ ) {
    $type = "script";
  } elsif( /^Makefile$/ or /^README$/ or /^TODO/ ) {
    $type = "text";
  } elsif( /^Rules.make$/ or /^Include.make$/ ) {
    $type = "text";
  } elsif( /\.txt$/ ) {
    $type = "text";
  } elsif( /\.str$/ ) {
    $type = "text";
  } else {
    my $output = qx/file $file/;
    $type = parse_file_output($output);
    #print STDERR "$type $output";
    #$type = "";
  }

  return $type;
}

sub parse_file_output
{
  my $output = shift;

  $output =~ s/^.+:\s*//;

  if( $output =~ /^Perl script/ ) {
    return "script";
  } elsif( $output =~ /^a.*perl.*script/ ) {
    return "script";
  } elsif( $output =~ /^a.*wish.*script/ ) {
    return "script";
  } elsif( $output =~ /^POSIX shell script/ ) {
    return "script";
  } elsif( $output =~ /^Bourne-Again shell script/ ) {
    return "script";
  } else {
    return "";
  }
}
