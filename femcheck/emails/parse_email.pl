#!/usr/bin/perl -w
#
# reads file addresses.txt, parses it and writes it to stdout
#
#--------------------------------------------------------

use strict;
use warnings;

%::alias = ();
my @addresses = ();
my $addresses;

read_alias();

$addresses = get_addresses("shyfem_g");
push(@addresses,@$addresses);
$addresses = get_addresses("shyfem_d");
push(@addresses,@$addresses);
$addresses = get_addresses("shyfem_u");
push(@addresses,@$addresses);

foreach my $address (@addresses) {
  print "$address\n";
}

#--------------------------------------------------------

sub read_alias {

  while(<>) {

    chomp;

    next if /^\s*$/;
    next if /^#/;

    my @f = split;

    my $alias = shift(@f);
    die "*** error: no alias... parsing error\n" unless $alias eq "alias";

    my $name = shift(@f);

    my $value = join(" ",@f);
    $value =~ s/,/ /g;
    $value =~ s/\s+/ /g;

    if( $value =~ /</ ) {
      $::alias{$name} = $value;
    } else {
      my @g = split(/\s+/,$value);
      $::alias{$name} = \@g;
    }
  } 

  my $error = 0;
  foreach my $key (keys %::alias) {
    my $value = $::alias{$key};
    if ( not defined $value ) {
      print STDERR "no such alias: $key\n";
      $error++;
    }
  }

  die "*** error: some alias could not be found\n" if $error;
}

#----------------------------------------------------

sub get_addresses {

  my $name = shift;

  my @array = ();
  my @addresses = ();
  my $error = 0;

  push(@array,$name);

  while( my $next = shift(@array) ) {
    my $value = $::alias{$next};
    if ( not defined $value ) {
      print STDERR "no such alias: $next\n";
      $error++;
    } elsif (ref $value eq 'ARRAY') {
      push(@array,@$value);
    } else {
      push(@addresses,$value);
    }
  }

  die "*** error: some alias not found\n" if $error;

  return \@addresses;
}

#----------------------------------------------------

