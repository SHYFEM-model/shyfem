#!/usr/bin/perl -s
#
# creates dependencies for fortran source files
#
# no recursive includes allowed -> changed: one level of recursion is allowed
#
# switch : -make

if( $make ) {
  #print "writing to makefile...\n"
}

while(<>) {

  if( /^\s+include\s*[\'\"]\s*(\S+)\s*[\'\"]\s*$/i ) {
    $name = $1;

    $inputfile = $ARGV;
    if( $inputfile =~ /^(\w+)\.(\w+)$/ ) {
      $file = $1;
      $ext = $2;
    } else {
      die "Cannot parse file name: $inputfile\n";
    }

    #print STDERR "$inputfile:  $file  $ext  $name\n";

    if( $ext eq "h" ) {
      $deph{$inputfile} .= "$name ";
    } else {
      $dep{$file} .= "$name ";
    }
  }

}

&open_makefile;

foreach $file (keys %dep) {

  $depfiles = $dep{$file};
  $uniq = &unique( $depfiles );
  $line = subst_h_files($uniq);
  
  #print "$file.o: $depfiles\n";
  #print "$file.o: $line\n";

  &write_makefile( "$file.o: $line" );

}

&close_makefile;

#################################################################

sub subst_h_files {

  my $line = shift;

  foreach $hfile (keys %deph) {
    if( $line =~ / $hfile / ) {
      $line .= $deph{$hfile};
      #print STDERR "** subst_h_files: $line\n";
    }
  }

  return $line;

}

#################################################################

sub unique {

  my $f = shift;
  my $old = "";
  my @new = ();

  my @f = split(/\s+/,$f);

  @f = sort( @f );

  foreach $f (@f) {
    if( $f ne $old ) {
      push(@new,$f);
      $old = $f;
    }
  }

  my $line = join(" ",@new);
  $line = " $line ";
  return $line;
}

#################################################################

sub open_makefile {

  if( $make ) {
    if( -f "makefile" ) {
      $makefile = "makefile";
    } elsif( -f "Makefile" ) {
      $makefile = "Makefile";
    } else {
      $makefile = "";
    }
  } else {
    $makefile = "";
  }

  if( $makefile ) {
    #print STDERR "opening $makefile\n";
    open(MAKE,">>$makefile");
    print MAKE "\n";
  }
}

sub close_makefile {

  if( $makefile ) {
    print MAKE "\n";
    #print STDERR "closing $makefile\n";
    close(MAKE);
  }
}

sub write_makefile {

  my $dep = shift;

  if( $makefile ) {
    print MAKE "$dep\n";
  } else {
    print "$dep\n";
  }
}

#################################################################

