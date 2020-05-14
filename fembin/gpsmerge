#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# merges ps files

$page = 0;
$totpages = 0;

while( $file = shift ) {

  $pages = &read_file( $file );
  $totpages += $pages;

  &handle_header;
  &handle_body;
  &handle_trailer;

}

&print_trailer;

######################################################

sub read_file {

  my $file = shift;
  my @content = ();

  open(INPUT,"$file");
  @content = <INPUT>;
  close(INPUT);

  @fileheader = ();
  @filebody = ();
  @filetrailer = ();
  my $where = 1;
  my $pages = 0;

  foreach (@content) {

    if( /^%%Page:/ ) { $where = 2; $pages++; }
    if( /^%%Trailer/ ) { $where = 3; }

    if( $where == 1 ) {
	push(@fileheader,$_);
    } elsif( $where == 2 ) {
	push(@filebody,$_);
    } else {
	push(@filetrailer,$_);
    }
  }

# if nothing in body, no DSC comments found -> everything is body

  my $n = @filebody;
  if( $n == 0 ) {
    die "Cannot merge PS file without DSC comments ($file)\n";
    @filebody = @fileheader;
    @fileheader = ();
    $pages = -1;
  }

  return $pages;
}

######################################################

sub handle_header {

  my $n = @fileheader;
  my $i;
  my $line;

  for($i=0;$i<$n;$i++) {
    $line = $fileheader[$i];
    if( $line =~ /^%%Pages:/ ) {
      $fileheader[$i] = "%%Pages: (atend)\n";
    }
    if( $line =~ /^%%BoundingBox:/ ) {
      &treat_bounding_box($line);
      $fileheader[$i] = "%%BoundingBox: (atend)\n";
    }
  }
  
  if( $nheader == 0 ) {		#print header only first time
    foreach (@fileheader) {
      print;
    }
  }

  $nheader++;
}

######################################################

sub handle_body {

  my $n = @filebody;
  my $i;
  my $line;

  for($i=0;$i<$n;$i++) {
    $line = $filebody[$i];
    if( $line =~ /^%%Page:/ ) {
      $page++;
      $filebody[$i] = "%%Page: $page $page\n";
    }
  }
  
  foreach (@filebody) {
    print;
  }
}

######################################################

sub handle_trailer {

  foreach (@filetrailer) {
    if( /^%%(\w+):\s*(.*)/ ) {
      my $key = $1;
      my $value = $2;
      $trailer{$key} = $value unless $trailer{$key};
      if( $key eq "BoundingBox" ) {
	&treat_bounding_box($_);
      }
    } elsif( /^%%(\w+)/ ) {
      my $key = $1;
      $trailer{$key} = "-" unless $trailer{$key};
    } elsif( $ntrailer = 0 ) {
      push(@trailertext,$_);
    }
  }

 $ntrailer++;
}

######################################################

sub print_trailer {

  $trailer{"Trailer"} = "";
  $trailer{"BoundingBox"} = "";
  $trailer{"Pages"} = "";
  $trailer{"EOF"} = "";

  print "%%Trailer\n";
  print "%%Pages: $page\n";
  if( $boundingbox ) {
    print "%%BoundingBox: $btot_llx $btot_lly $btot_urx $btot_ury\n";
  }
  foreach (keys %trailer) {
    $val = $trailer{$_};
    if( $val eq "" ) {
      ;
    } elsif( $val eq "-" ) {
      print "%%$_\n";
    } else {
      print "%%$_: $val\n";
    }
  }
  print "%%EOF\n";
}

######################################################

sub treat_bounding_box {

  my $line = shift;

  if( $line =~ /^%%BoundingBox:\s+(.*)/ ) {
    $bb = $1;
    #print "%%PageBoundingBox: $_\n";
    unless( $bb =~ /\(atend\)/ ) {
	@bb = split("\s+",$bb);;
	$btot_llx = $bb[0] if( $btot_llx eq "" || $btot_llx > $bb[0] );
	$btot_lly = $bb[1] if( $btot_lly eq "" || $btot_lly > $bb[1] );
	$btot_urx = $bb[2] if( $btot_urx eq "" || $btot_urx < $bb[2] );
	$btot_ury = $bb[3] if( $btot_ury eq "" || $btot_ury < $bb[3] );
    }
  }

  $boundingbox++;

  return $bb;

}





