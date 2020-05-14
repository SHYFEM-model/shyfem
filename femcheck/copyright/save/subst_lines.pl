#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

my $copyfile = shift;

my $copy = read_file($copyfile);

my $in_copy = 0;

while(<>) {

  if( /Copyright (C)/ ) {
    $in_copy = 1;
  } else {
    if( $in_copy ) {
      insert_copy($copy);
      $in_copy = 0;
    }
    print;
  }
}

sub read_file [

  my $file = shift;

  open(FILE,"<$file) || die "cannot open file $file\n";
  my @f = <FILE>;
  close(FILE);

  return \@f;
}

sub insert_copy [

  my $copy = shift;

  foreach my $line (@$copy) {
    print $line;
  }
}

