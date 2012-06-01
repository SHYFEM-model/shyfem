#!/usr/bin/perl -w
#
# grd utility routines
#
# example of usage:
#
# use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
#
# use csv;
# 
# my $file = $ARGV[0];
# my $csv = new csv;
# $csv->readfile($file);
#
# or
#
# my $csv = new csv($file);
#
# if split on whitespace, use 
#
# my $csv = new csv($file,"\\s+");
#
# or 
#
# my $csv = new csv($file,'\s+');
# 
# columns are numbered from 1 to n
# in col[0] is the time column (maybe scaled)
# readfile resets all variables for the file (number of cols, etc..)
#
##############################################################
#
# 1.1	30.05.2012	ggu	some new functionality
# 1.0	01.12.2005	ggu	bug fixed in writefile -> novalue for undefined
#
##############################################################

use strict;

package csv;

##############################################################

sub new
{
    my ($pck,$file,$separator) = @_;

    my $self =	{
	    		 file		=>	undef
			,lines		=>	undef
			,rows		=>	undef
			,cols		=>	undef
			,names		=>	undef
			,separator	=>	";"
			,novalue	=>	"-999"
			,extension	=>	"dat"
			,prefix		=>	""
			,nlines		=>	0
			,ncols		=>	0
			,write_time	=>	1
			,verbose	=>	1
			,version	=>	1.1
		};

    bless $self;

    $self->readfile($file,$separator);

    return $self;
}

####################################################################

sub readfile {

  my ($self,$file,$separator) = @_;

  return unless $file;
  open(FILE,"$file") || die "Cannot open file: $file\n";

  $self->{file} = $file;
  $self->{separator} = $separator if $separator;
  $separator = $self->{separator};

  #print STDERR "separator: $separator\n";

  my @rows = ();
  my @lines = ();
  my $ncols = 0;

  while( my $line = <FILE> ) {

	$line =~ s/\r*\n//;
	$line =~ s/^\s+//;			#get rid of leading blanks
	my @f = split(/$separator/,$line);

	my $n = @f;
	if( $n > $ncols ) {
	  $ncols = $n;
	}

	push(@rows,\@f);
	push(@lines,$line);

	#print STDERR "$line\n";
  }

  close(FILE);

  my $nlines = @lines;
  $self->{nlines} = $nlines;
  $self->{ncols} = $ncols;
  $self->{rows} = \@rows;
  $self->{lines} = \@lines;

  #print STDERR "$nlines $ncols\n";

  my @cols = ();
  for(my $i=1;$i<=$ncols;$i++) {
    my @col = ();
    $cols[$i] = \@col;
  }

  #print STDERR "$nlines $ncols\n";

  foreach my $row (@rows) {
    for(my $i=1;$i<=$ncols;$i++) {
      my $col = $cols[$i];
      my $val = $row->[$i-1];
      if( not defined($val) or $val eq "" ) {
        $val = $self->{novalue};
      }
      push(@$col,$val);
    }
  }

  #print STDERR "$nlines $ncols\n";

  $self->{cols} = \@cols;
  $self->set_time_column(0);

  if( $self->{verbose} ) {
    print STDERR "reading file: $file ... records: $nlines  columns: $ncols \n";
  }
}

#-------------------------------------------------------

sub writefile {

  my ($self,$file,@icols) = @_;

  if( $file ) {
    my $ext = $self->{extension};
    my $pre = $self->{prefix};
    $file = "$pre$file.$ext";
    open(FILE,">$file") || die "Cannot open file: $file\n";
  } else {
    die "No file name given\n";
  }

  my $debug = 0;

  my $write_time = $self->{write_time};
  my $novalue = $self->{novalue};
  my $nlines = $self->{nlines};
  my $ncols = $self->{ncols};
  my $cols = $self->{cols};
  my $time = $cols->[0];

  @icols = (2..$ncols);

  for(my $i=0;$i<$nlines;$i++) {
      print FILE "$time->[$i] " if $write_time;
      foreach my $ic (@icols) {
        my $col = $cols->[$ic];
        my $val = $col->[$i];
	unless( defined $val ) {
	  print STDERR "$file $time->[$i] $i $ic - no value\n" if $debug;
	  $val = $novalue;
	}
        print FILE "$val ";
      }
      print FILE "\n";
  }

  close(FILE);
}

sub writefiles {

  my ($self,$basename) = @_;

  my $ncols = $self->{ncols};
  my $extension = $self->{extension};
  my $file;

  for(my $i=1;$i<=$ncols;$i++) {
    if( $basename ) {
      $file = "${basename}_$i";
    } else {
      $file = $self->get_filename($i);
    }
    print "$i $file\n";
    $self->writefile($file,$i);
  }
}
  
sub info {

  my ($self) = @_;

  my $names = $self->{names};

  print "file name: $self->{file}\n";
  print "separator: $self->{separator}\n";
  print "extension: $self->{extension}\n";
  print "lines:     $self->{nlines}\n";
  print "columns:   $self->{ncols}\n";

  print "names of columns:\n";
  my $sep = "";
  foreach my $item (@$names) {
    print "$sep$item";
    $sep = "|";
  }
  print "\n";
}

#-------------------------------------------------------

sub get_name {

  my ($self,$icol) = @_;

  my $names = $self->{names};

  return $names->[$icol];
}

sub get_filename {

  my ($self,$icol) = @_;

  my $name = $self->get_name($icol);

  $name =~ s/\(.*$//;	#eliminate pars
  $name =~ s/^\s+//;	#trim white space
  $name =~ s/\s+$//;	#trim white space
  $name =~ s/\s+/_/g;	#substitute space

  $name = "column_$icol" unless $name;

  return $name;
}

sub set_header {

# $head		total number of header rows
# $name		row which contains names of columns

  my ($self,$head,$name) = @_;

  my $ncols = $self->{ncols};
  my $rows = $self->{rows};
  my $lines = $self->{lines};
  $self->{nlines} -= $head;
  my $header_names = "";

  for(my $i=1;$i<=$head;$i++) {
    my $save = shift(@$rows);
    $header_names = $save if $i == $name;
  }

  unless ($header_names) {	# no header names are given
    my @f = ();
    for(my $i=1;$i<=$ncols;$i++) {
      $name = "column_$i";
      push(@f,$name);
    }
    $header_names = \@f;
  }

  unshift(@$header_names,"time");
  $self->{names} = $header_names;

  my $cols = $self->{cols};

  foreach my $col (@$cols) {
    for(my $i=1;$i<=$head;$i++) {
      shift(@$col);
    }
  }

  $self->set_time_column(0);
}

sub extend_data {

  my ($self,$infront,$cycle,$dt) = @_;

  my $cols = $self->{cols};

  my $isrc;
  if( ( $infront and $cycle ) or ( not $infront and not $cycle ) ) { 
    $isrc = -1;		#last
  } else {
    $isrc = 0;		#first
  }

  #$self->print_col(0);

  foreach my $col (@$cols) {
    my $value = $$col[$isrc];
    if( $infront ) {
      unshift(@$col,$value);
    } else {
      push(@$col,$value);
    }
  }

  #$self->print_col(0);

  $self->{nlines}++;
  #$self->{nlines} += 1;

  # adjust time column

  my $col = $$cols[0];
  if( $infront ) {
      $$col[0] = $$col[1] - $dt;
  } else {
      $$col[-1] = $$col[-2] + $dt;
  }

}

sub print_col {

  my ($self,$icol) = @_;

  my $col = $self->{cols}->[$icol];

  print "Printing column $icol\n";
  foreach my $item (@$col) {
    print "$item\n";
  }
}

sub scale_col {

  my ($self,$icol,$fact,$v0) = @_;

  my $col = $self->{cols}->[$icol];
  $v0 = 0 unless $v0;

  foreach my $item (@$col) {
    $item = $v0 + $fact * $item;
  }
}

sub get_col {

  my ($self,$col) = @_;

  my $cols = $self->{cols};

  return $cols->[$col];
}

sub copy_col {

  my ($self,$col) = @_;

  my $cols = $self->{cols};
  my $column = $cols->[$col];
  my @f = ();

  foreach my $item (@$column) {
    push(@f,$item);
  }

  return \@f;
}

sub set_time_column {

  my ($self,$col,$fact,$t0) = @_;

  my $time;
  if( $col > 0 ) {
    $time = $self->copy_col($col);
  } else {
    $time = $self->make_range( $self->{nlines} );
  }
  $self->{cols}->[0] = $time;

  return if not $fact and not $t0;

  $fact = 1 unless $fact;
  $t0 = 0 unless $t0;
  $self->scale_col(0,$fact,$t0);
}

sub scale_time_column {

  my ($self,$fact,$t0) = @_;

  $self->scale_col(0,$fact,$t0);
}

sub make_range {

  my ($self,$max) = @_;

  my @array = ();
  for(my $i=1;$i<=$max;$i++) {
    push(@array,$i);
  }
  return \@array;
}
  
sub set_separator {

  my ($self,$sep) = @_;

  $self->{separator} = $sep;
}

sub set_novalue {

  my ($self,$novalue) = @_;

  $self->{novalue} = $novalue;
}

sub set_verbose {

  my ($self,$verbose) = @_;

  $self->{verbose} = $verbose;
}

sub set_extension {

  my ($self,$ext) = @_;

  $self->{extension} = $ext;
}

sub set_prefix {

  my ($self,$pre) = @_;

  $self->{prefix} = $pre;
}

sub no_time_column {	#call if no time column present in file

  my ($self) = @_;

  $self->{write_time} = 0;
}

#-------------------------------------------------------

###################################

# example for using outhandle -> delete afterwards

sub write_node
{
    my ($self,$item) = @_;

    my $oldhandle = undef;

    my $handle = $self->{outhandle};
    if( $handle ) {
      $oldhandle = select $handle;
    }

    print "1 $item->{number} $item->{type}";
    print " $item->{x} $item->{y}";
    print " $item->{h}\n" if $item->{h};
    print "\n";

    select $oldhandle if $oldhandle;
}

###################################

###################################
#&test_grd(@ARGV);
###################################
1;
###################################

