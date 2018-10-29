#!/usr/bin/perl -w
#
#-----------------------------------

use strict;

my $copyright = "/home/georg/shyfem/femcheck/copyright/copyright_notice.txt";

my $lang = shift;
$lang = "" unless $lang;

if( $lang eq "f" ) {
  ;
} elsif( $lang eq "c" ) {
  ;
} elsif( $lang eq "t" ) {
  ;
} else {
  print "Usage: include_copyright.pl {f|c|t} file\n";
  die "unknown language: $lang\n";
}

my $copy = read_copyright();
$copy = adjust_copyright($copy,$lang);
print_copyright($copy);

print_file($lang);

#-----------------------------------

sub adjust_copyright
{
  my ($copy,$lang) = @_;

  my ($char,$bline,$eline,$line);

  if( $lang eq "f" ) {
    $char = "!";
    my $stars = "-" x 74;
    $bline = "!" . $stars . "\n";
    $eline = $bline;
  } elsif( $lang eq "c" ) {
    $char = " *";
    my $stars = "*" x 72;
    $bline = "/" . $stars . "\\\n";
    $eline = "\\" . $stars . "/\n";
  } elsif( $lang eq "t" ) {
    $char = "%";
    my $stars = "-" x 72;
    $bline = "%" . $stars . "\n";
    $eline = $bline;
  } else {
    die "unknown language: $lang\n";
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
  open(COPY,"<$copyright") || die "cannot read file: $copyright\n";
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
  my $lang = shift;

  my $debug = 0;	# set to > 0 for debug
  my $header = 1;	# we are still in header
  my $copy = 0;		# we are in old copyright

  while(<>) {
    if( $header ) {
      if( $lang eq "f" ) {
        next if /^[cC!]\s*$/;
        next if /^\s*$/;
        next if /^[cC!]\s+\$Id:/;
      } elsif( $lang eq "c" ) {
	next if /^\/\* \$Id:/;
        $copy = 1 if /Copyright \(c\)/;
        if( /Revision / ) {
          $copy = 0;
	  $header = 0;
        }
	next if $copy;
	print;
	next;
      } elsif( $lang eq "t" ) {
        next if /^\s*$/;
        next if /^%\s+\$Id:/;
      } else {
	die "unknown language: $lang\n";
      }
    }
    last if( --$debug == 0 );	# only executed if debug was > 0
    print;
    $header = 0;
  }
}

