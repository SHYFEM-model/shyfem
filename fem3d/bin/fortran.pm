#!/usr/bin/perl -ws
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# parses fortran routines
#
# version 1.0			written from scratch
# version 1.1			re-written, but cannot do structure
# version 1.2	02.12.2014	adapted for computing structure
# version 1.3	08.12.2014	fully functional
#
#---------------------------------------------------------
#
# example of usage:
#
# use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
#
# use fortran;
# 
# my $fortran = new fortran;
#
# $fortran->set_ignore_files(\@ignore_files);
# $fortran->parse_files($quiet);
# $fortran->show_functions();
# $fortran->check_uniqueness();
# $fortran->parse_routines();
#
# still to do:
#
#	allow for statement function in midst specification statements
#	allow of comments inside continuation line
#	allow for ! inside string
#	should reject real*8
#	should reject assigned goto
#	equivalence:, subrgf.f, sublnka.f, debug.f
#
#########################################################

use strict;

package fortran;

#########################################################

sub new
{
	my $self;

	$self = 	{
				 defined	    =>	{}
				,functions	    =>	{}
				,called		    =>	{}

				,calling	    =>	{}
				,calledby	    =>	{}
				,defined_files	    =>	{}

				,all_routines	    =>	{}
				,routines	    =>	{}
				,files		    =>	{}

				,ignore_files	    =>	{}
				,ignore_routines    =>	{}

				,debug		    =>	0
				,parse_later	    =>	0
				,upper		    =>	0
				,no_main	    =>	0

				,comnum		    =>	0
				,implicit	    =>	0
				,contains	    =>	0

				,file		    =>	""
				,nlines		    =>	0
				,ignore		    =>	0
				,has_program	    =>	0
				,code		    =>	[]
				,space		    =>	""

				,calling_names	    =>	{}
				,calling_times	    =>	{}
				,calling_stack	    =>	[]

				,maxl		    =>	0	#actual
				,maxlevel	    =>	0	#permitted
				,maxlevel_reached   =>	0	#how often over
				,level		    =>	0	#actual level
			};

	bless $self;
	return $self;
}

#%defined = ();		# how often is this routine defined
#%functions = ();	# how often is this function defined
#%called = ();		# how often is this routine called

#%calling = ();		# routine calls these routines
#%calledby = ();	# routine is called by these routines
#%defined_files = ();	# files where routine has been defined

#%routines = ();	# complete info on routines
#%files = ();		# complete info on file (fitem)

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub parse_files {

 my ($self,$quiet) = @_;

  foreach my $file (@ARGV) {
    my $ignore = $self->{ignore_files}->{$file};
    if( $ignore ) {
      print STDERR "ignoring $file\n" unless $quiet;
    } else {
      print STDERR "reading $file\n" unless $quiet;
      $self->parse_file($file,$quiet);
      $self->parse_routine_by_file($file) unless $self->{parse_later};
    }
  }
}

sub parse_file {

  my ($self,$file,$quiet) = @_;

  my $ignore = 0;
  my $in_interface = 0;
  my ($comm,$next,$name,$type,$code);

  open(FILE,"<$file") || die "Cannot open file: $file\n";

  $self->open_file($file);

  while(1) {
    ($comm,$next) = $self->read_comment();
    $name = $self->make_section_name("comment");
    $self->insert_routine($name,'comment',$comm);
    last unless defined $next;
    ($name,$type,$code) = $self->read_routine($next);
    $self->insert_routine($name,$type,$code);
  }

  $self->close_file();
  close(FILE);
}

sub read_comment {

  my ($self) = @_;

  my @comm = ();

  while(<FILE>) {
    chomp;
    $self->{nlines}++;
    #print STDERR "in comment: $_\n";
    last if( not /^\s*$/ and not /^[cC\*\!]/ and not /^\s*\!/ );
    push(@comm,$_);
  }

  return (\@comm,$_);
}

sub make_section_name {

  my ($self,$name) = @_;

  $self->{comnum}++;
  return $name . "_" . $self->{comnum};
}

sub insert_routine {

  my ($self,$name,$type,$code) = @_;

  my $n = @$code;
  return unless $n;

  my $file = $self->{file};

  my $ritem = {};
  $ritem->{name} = $name;
  $ritem->{type} = $type;
  $ritem->{file} = $file;
  $ritem->{ignore} = $self->{ignore};
  $ritem->{contain} = $self->{contain};
  $ritem->{code} = $code;

  my $fitem = $self->{files}->{$file};
  my $seq = $fitem->{sequence};
  push(@$seq,$ritem);

  return if $type eq "comment";
  return if $type eq "end";
  $self->{functions}->{$name}++ if $type eq "function";
  $self->{all_routines}->{"$name-$file"} = $ritem;

  if( $self->{no_main} == 0 or $self->{has_program} == 0 ) {
    $self->{routines}->{$name} = $ritem;
    $self->{defined}->{$name}++;
    $self->{defined_files}->{$name} .= "$file ";
  }
}

sub read_routine {

  my ($self,$next) = @_;

  my $in_interface = 0;
  my $orig = "";

  my @code = ();
  push(@code,$next);

  my ($name,$type) = $self->parse_initial($next);

  if( $type eq 'end' ) {
    $name = "" unless $name;
    my $cont = $self->{contains};
    unless( $cont ) {
      die "end found but no contains ... aborting\n";
    }
    print STDERR "closing contains statement of $cont\n" if $self->{debug};
    $self->{contains}--;
    $name = "end" . "_" . $cont;
    $name = $self->make_section_name($name);
    return ($name,$type,\@code);
  }

  #print STDERR "reading routine: $name $type $next\n";

  while (<FILE>) {

    chomp;
    $self->{nlines}++;
    $orig = $_;
    push(@code,$_);

    s/^\s*\d+\s+/ /;	# get rid of label number
    s/\s+//g;		# no space
    #s/\!.*//g;		# line comment
    $_ = $self->clean_comment($_) if /\!/;

    if( /^INTERFACE$/i ) { 
      $in_interface = 1;
    } elsif( /^ENDINTERFACE$/i ) { 
      $in_interface = 0;
    }
    next if $in_interface;
    
    if( /^CONTAINS$/i ) { 
      $self->{contains}++;
      last;
    }

    last if( /^END$type$name/i );
    last if( /^END$type/i );
    last if( /^END$/i );
  }

  if ( not /^END/i and not /^CONTAINS/i ) {
    die "*** cannot find end statement of $type $name: $orig\n";
  }

  return ($name,$type,\@code);
}

sub clean_comment {

  my ($self,$string) = @_;

  my $orig = $_;
  my $ap = 0;
  my @new = ();

  foreach my $c (split //, $string) {
    $ap++ if $c eq "'";
    if( $c eq '!' ) {
      last if $ap%2 == 0;
    }
    push(@new,$c);
  }

  $string = join("",@new);

  #print STDERR "clean_comment: |$orig|$_|\n";

  return $string;
}

sub parse_initial {

  my ($self,$next) = @_;

  $_ = $next;
  s/^\s*\d+\s+/ /;	# get rid of label number
  s/\s+//g;		# no space
  #s/\!.*//g;		# line comment
  $_ = $self->clean_comment($_) if /\!/;

  if( /^PROGRAM(\w+)/i ) { 
    die "*** more than one main program: $next\n" if $self->{has_program};
    $self->{has_program}=1;
    return ($1,'program');
  } elsif( /^SUBROUTINE(\w+)/i ) { 
    return ($1,'subroutine');
  } elsif( /^FUNCTION(\w+)/i ) { 
    return ($1,'function');
  } elsif( /^BLOCKDATA(\w+)/i ) { 
    return ($1,'blockdata');
  } elsif( /^(REAL|INTEGER|DOUBLEPRECISION|LOGICAL)FUNCTION(\w+)/i ) { 
    return ($2,'function');
  } elsif( /^MODULE(\w+)/i ) { 
    return ($1,'module');
  } elsif( /^END(\w*)/i ) { 
    return ($1,'end');
  } else {
    print STDERR "*** main definition without program statement: $next\n";
    $self->{has_program}=1;
    return ('main','program');
  }
}

sub open_file {

  my ($self,$newfile) = @_;

  my %fitem = ();
  $self->{files}->{$newfile} = \%fitem;

  $fitem{filename} = $newfile;
  $fitem{changed} = 0;
  $fitem{nlines} = 0;
  $fitem{has_program} = 0;
  $fitem{sequence} = [];;

  $self->{file} = $newfile;
  $self->{nlines} = 0;
  $self->{has_program} = 0;	# is 1 if file contains program statement
}

sub close_file {

  my ($self) = @_;

  return unless $self->{file};

  my $fitem = $self->{files}->{$self->{file}};
  $fitem->{nlines} = $self->{nlines};
  $fitem->{has_program} = $self->{has_program};
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub get_routine_list {

  my ($self) = @_;

  return sort keys %{$self->{routines}};
}

sub get_routine_item {

  my ($self,$routine) = @_;

  return $self->{routines}->{$routine};
}

sub get_file_list {

  my ($self) = @_;

  return sort keys %{$self->{files}};
}

sub get_file_item {

  my ($self,$filename) = @_;

  return $self->{files}->{$filename};
}

#----------------------------------------------------------------

sub show_files {

  my ($self) = @_;

  print STDERR "list of files:\n";
  foreach my $file (sort keys %{$self->{files}}) {
    my $fitem = $self->{files}->{$file};
    my $filename = $fitem->{filename};
    my $sequence = $fitem->{sequence};
    print "  $filename $file \n";
    foreach my $ritem (@$sequence) {
      my $name = $ritem->{name};
      my $type = $ritem->{type};
      print "    $name ($type)\n";
    }
  }
}

sub show_routines {

  my ($self) = @_;

  print STDERR "list of routines:\n";
  my @list = sort keys %{$self->{routines}};
  $self->show_procedures(\@list);
}

sub show_functions {

  my ($self) = @_;

  print STDERR "list of functions:\n";
  my @list = sort keys %{$self->{functions}};
  $self->show_procedures(\@list);
}

sub show_procedures {

  my ($self,$list) = @_;

  foreach my $name (@$list) {
    my $count = $self->{defined}->{$name};
    $count = 0 unless $count;
    my $item = $self->{routines}->{$name};
    my $file = $item->{file};
    print "  $count  $name  $file\n";
  }
}

sub check_uniqueness {

  my ($self,$more) = @_;

  foreach my $name (keys %{$self->{defined}}) {
    my $ritem = $self->{routines}->{$name};
    next if $ritem->{ignore};
    my $count = $self->{defined}->{$name};
    if( $count != 1 ) {
      print STDERR "multiple definition of routine: $count  $name\n";
      if( $more ) {
	my $line = $self->{defined_files}->{$name};
        print STDERR "  $line\n";
      }
    }
  }
}

sub show_content {

  my ($self) = @_;

  print STDERR "content of routines:\n";

  my @list = sort keys %{$self->{routines}};
  foreach my $rname (@list) {
    my $ritem = $self->{routines}->{$rname};
    my $name = $ritem->{name};
    my $file = $ritem->{file};
    my $implicit = $ritem->{implicit};
    print "  $name: $implicit ($file)\n";
    my $common = $ritem->{common};
    foreach $name (sort keys %$common) {
      my $list = $common->{$name};
      my $line = join(",",@$list);
      print "    common /$name/ $line\n";
    }
    my $save = $ritem->{save};
    if( $save ) {
      my $line = join(",",sort keys %$save);
      print "    save $line\n";
    }
    my $declaration = $ritem->{declaration};
    foreach my $type (sort keys %$declaration) {
      my $htype = $ritem->{declaration}->{$type};
      my $line = join(",",sort keys %$htype);
      print "    $type $line\n";
    }
  }
}

sub show_implicit {

  my ($self) = @_;

  return unless $self->{implicit};

  print STDERR "content of files:\n";
  foreach my $file (sort keys %{$self->{files}}) {
    my $fitem = $self->{files}->{$file};
    my $filename = $fitem->{filename};
    my $implicit = $fitem->{implicit};
    print "  $filename: $implicit \n" if $implicit;
  }
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub parse_routine_by_file {

  my ($self,$file) = @_;

  $::ntotlines = 0;
  my $fitem = $self->{files}->{$file};
  my $sequence = $fitem->{sequence};

  foreach my $ritem (@$sequence) {
    my $name = $ritem->{name};
    my $type = $ritem->{type};
    my $code = $ritem->{code};
    my $ncode = @$code;
    if( $type ne "comment" and $type ne "end" ) {
      $self->parse_routine($ritem);
    }
    $::ntotlines += $ncode;
  }
}

sub parse_routines {

  my ($self) = @_;

  foreach my $name (keys %{$self->{defined}}) {
    my $ritem = $self->{routines}->{$name};
    $self->parse_routine($ritem);
  }
}

sub parse_routine {

  my ($self,$ritem) = @_;

  my $name = $ritem->{name};
  my $type = $ritem->{type};
  my $file = $ritem->{file};
  my $code = $ritem->{code};

  my $in_executable = 0;

  #print STDERR "parsing $name ($type)\n";

  unless( $code ) {
    print STDERR "*** no code for routine $name\n";
    return;
  }

  my $cline = 0;
  my $nline = 0;
  my $n = @$code;

  for(my $i=0; $i<$n; $i++) {

    $cline++;
    $nline++;
    my $line = $code->[$i];
    $line = $self->clean_line($line);
    next unless $line;
    while( $i+1 < $n and $code->[$i+1] =~ /^     \S/ ) {	# continuation
      $i++; $nline++;
      my $aux = $code->[$i];
      $aux = $self->clean_line($aux);
      #print "conti found... $aux\n";
      $aux =~ s/^     \S/ /;
      $line .=  $aux;
    }

    next if( $cline == 1 );		# skip opening of routine

    print STDERR "parsing: $line\n" if $::bdebug;
    if( not $in_executable ) {
      print STDERR "checking spec\n" if $::bdebug;
      next if $self->is_specification($line,$ritem);
      print STDERR "not spec\n" if $::bdebug;
      $in_executable = 1;
    }
    if( $in_executable ) {
      print STDERR "checking exec\n" if $::bdebug;
      unless( $self->is_executable($line,$ritem) ) {
	my $nnline = $::ntotlines + $nline;
	if( $self->is_specification($line,$ritem) ) {
          #die "*** specification in executable ($name,$nline,$file): $line\n";
          die "*** specification in executable ($nnline): $line\n";
	} else {
          #die "*** unknown statement ($name,$nline,$file): $line\n";
          die "*** unknown statement ($nnline): $line\n";
        }
      }
    }

    if( $line =~ /^\s+CALL\s+(\w+)/i ) { 
      $self->handle_call($name,$1);
    } elsif( $line =~ /^\s+IF\s*\(.+\)\s*CALL\s+(\w+)/i ) { 
	$self->handle_call($name,$1); 
    }
    while( $line =~ /(\w+)\s*\(/ig ) { 
      my $func = $self->upperlow($1);
      if( $self->{functions}->{$func} ) { 
	#print STDERR "function call: $func in routine $name\n";
        $self->handle_call($name,$func);
      }
    }
  }

  if( $self->{implicit} and not $ritem->{implicit} ) {
    print STDERR "no implicit none in $name (file $file)\n";
    my $fitem = $self->{files}->{$file};
    $fitem->{implicit} .= " $name ";
  }
}

sub clean_line {

  my ($self,$line) = @_;

  $line =~ s/^\s*\d+\s+/ /;		# get rid of label number
  $line =~ s/^[cC!*].*$//;		# line comment
  #$line =~ s/\s*!.*$//;		# trailing comment
  $line = $self->clean_comment($line) if $line =~ /\!/;
  $line =~ s/\s+$//;			# trailing space
  $line =~ s/^\s+$//;			# empty line

  return $line;
}

sub is_specification {

  my ($self,$line,$ritem) = @_;

  my $what = "found";
  my $type;

  $_ = $line;
  s/\s+//g;			# eliminate all space

  if( $self->is_assignment($_,$ritem) ) {
    return 0;
  } elsif( /^implicitnone\b/i ) {
    $ritem->{implicit} = 1;
  } elsif( /^implicit\w+\b/i ) { ;
  } elsif( /^include\'(\S+)\'/i or /^include\"(\S+)\"/i ) {
    $ritem->{include}->{$1} = 1;
  } elsif( /^use(\w+)\b/i ) {
    $ritem->{use}->{$1} = 1;
  } elsif( $type = $self->parse_type_declaration($_,$ritem) ) {
    $what = "declaration-$type";
  } elsif( /^common\//i ) { 
    $self->parse_common($_,$ritem); $what = "common";
  } elsif( /^equivalence\(/i ) { ;
  } elsif( /^parameter\(/i ) { ;
  } elsif( /^optional\w+/i ) { ;
  } elsif( /^save\w+\b/i or /^save\b/i ) { 
    $self->parse_save($_,$ritem); $what = "save";
  } elsif( /^data\w+\b/i ) { ;
  } elsif( /^dimension\w+\b/i ) { ;
  } elsif( /^external\w+\b/i or /^intrinsic\w+\b/i ) { ;
  } elsif( /^type\w+\b/i or /^type\b.*::\w+/i) { ;
  } elsif( /^endtype\w*\b/i ) { ;
  } elsif( /^interface\w*$/i or /^endinterface$/i ) { ;
  } elsif( /^subroutine\w+\b/i or /^endsubroutine\w*\b/i ) { ;
  } elsif( /^contains\b/i ) { ;
  } elsif( /^moduleprocedure\w+\b\b/i ) { ;
  } else {
    return "";
  }
  print STDERR "is spec: $_\n" if $::bdebug;
  return $what;
}

sub is_executable {

  my ($self,$line,$ritem) = @_;

  $_ = $line;
  s/\s+//g;			# eliminate all space

  if( $self->is_assignment($_,$ritem) ) { ;
  } elsif( /^call\w+\(/i or /^call\w+$/i ) { ;
  } elsif( /^if\(/i or /^elseif\(/i or /^else\b/i ) { ;
  } elsif( /^do\d*,?\w+=/i ) { ;
  } elsif( /^dowhile\(/i or /^do$/i ) { ;
  } elsif( /^continue\b/i ) { ;
  } elsif( /^write\(/i or /^read\(/i or /^print\*,/i ) { ;
  } elsif( /^open\(/i or /^close\(/i or /inquire\(/i ) { ;
  } elsif( /^rewind\(/i or /^rewind\w+$/i ) { ;
  } elsif( /^backspace\(/i or /^backspace\w+$/i ) { ;
  } elsif( /^endfile\(/i or /^endfile\w+$/i ) { ;
  } elsif( /^format\(/i or /^flush\(/ ) { ;
  } elsif( /^namelist\//i ) { ;
  } elsif( /^end\b/i or /^enddo\b/i or /^endif\b/i ) { ;
  } elsif( /^endsubroutine\w*\b/i or /^endfunction\w*\b/i ) { ;
  } elsif( /^endprogram\w*\b/i or /^endmodule\w*\b/i ) { ;
  } elsif( /^endtype\w*\b/i ) { ;
  } elsif( /^goto\d+\b/i or /^goto\(.+\),?\w+/i ) { ;
  } elsif( /^exit$/i or /^cycle$/i ) { ;
  } elsif( /^return\b/i or /^stop\b/i or /^stop\d+\b/i ) { ;
  } elsif( /^allocate\(/i or /^deallocate\(/i ) { ;
  } elsif( /^selectcase\(/i or /^endselect$/i ) { ;
  } elsif( /^case\(/i or /^casedefault$/i ) { ;
  } else {
    return 0;
  }
  print STDERR "is exec: $_\n" if $::bdebug;
  return 1;
}

sub is_assignment {

  my ($self,$line,$ritem) = @_;

  $line = elim_string($line);
  $line = elim_pars($line);
  print "assignment: $line\n" if $::bdebug;

  if( $line =~ /^(.*?)=(.*)$/ ) {
    my $ass = $1;
    my $val = $2;
    unless( defined $val ) { return 0; }
    if( $val =~ /,/ ) { return 0; }
    unless( defined $ass ) { return 0; }
    if( $ass =~ /^\w+$/ or $ass =~ /^\w+%\w+$/ ) { return 1; }
    #unless( $ass =~ /^\w+$/ ) { return 0; }
    return 0;
  } else {
    return 0;
  }
  return 1;
}

sub parse_type_declaration {

  my ($self,$line,$ritem) = @_;

  my @types = ('integer','real','logical','doubleprecision'
		,'complex','character','');

  my $orig = $line;
  my $print = 0;
  my ($type,$param,$dim,$list);
  $type = "";
  my $found = 0;

  foreach my $t (@types) {
    #print "trying $t: $line\n";
    last unless $t;
    if( $line =~ /^$t(.*)::(.+)$/i ) {			# real, save ::
      $type = $t; $param = $1; $list = $2; $found = 1; last;
    } elsif( $line =~ /^$t(\*\d+)(.+)$/i ) {		# real*8, character*80
      $type = $t; $dim = $1; $list = $2; $found = 2; last;
    } elsif( $line =~ /^$t(\*\(\*\))(.+)$/i ) {		# character*(*)
      $type = $t; $dim = $1; $list = $2; $found = 3; last;
#    } elsif( $line =~ /^$t(\*\(\d+\))(.+)$/i ) {	# character*(80)
#      $type = $t; $dim = $1; $list = $2; $found = "2a"; last;
#    } elsif( $line =~ /^$t(\(\*\))(.+)$/i ) {		# character(*)
#      $type = $t; $dim = $1; $list = $2; $found = "3a"; last;
    } elsif( $line =~ /^$t(\w.*)/i ) {			# real val
      $type = $t; $list = $1; $found = 4; last;
    }
  }
  #print "type found ($found): $type\n";

  return $type unless $type;
  my $alist = split_list($list);
  foreach my $name (@$alist) {
    $ritem->{declaration}->{$type}->{$name}++;
  }
  return $type;
}

sub parse_save {

  my ($self,$line,$ritem) = @_;

  my $orig = $line;
  $line =~ s/^save//;

  my $alist = split_list($line);

  foreach my $name (@$alist) {
    $ritem->{save}->{$name}++;
  }
}

sub parse_common {

  my ($self,$line,$ritem) = @_;

  my $orig = $line;
  $line =~ s/^common//;
  my ($name,$list);

  while(1) {
    ($name,$list,$line) = $self->next_common($line);
    last unless $name;
    $ritem->{common}->{$name} = split_list($list);
  }

  if( $line ) {
    die "*** cannot parse common block: $orig ($line)\n";   
  }
}

sub next_common {

  my ($self,$line) = @_;

  if( $line =~ /^\/(\w+)\/([^\/]+)(.*)/ ) {
    my $name = $1;
    my $list = $2;
    $line = $3;
    $list =~ s/,$//;	# pop trailing comma
    return ($name,$list,$line);
  }
  return;
}

sub split_var {

  my ($self,$line) = @_;

  return split_list($line);
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub split_list {

  my $list = shift;

  return [] unless $list;

  my @alist = ();

  my $m = 0;
  my $new = "";
  my @f = split(//,$list);
  foreach my $c (@f) {
    if( $c eq "," and $m == 0 ) {
      push(@alist,$new);
      $new = "";
    } else {
      $new .= "$c";
    }
    $m++ if $c eq '(';
    $m-- if $c eq ')';
    die "*** cannot parse parenthesis: $list\n" if $m < 0;
  }
  die "*** cannot parse parenthesis: $list\n" if $m > 0;
  push(@alist,$new) if $new;

  return \@alist;
}

sub elim_string {

  my $token = shift;

  return $token unless $token =~ /[\'\"]/;

  my $close = "";
  my $new = "";
  my @f = split(//,$token);
  foreach my $c (@f) {
    if( $c eq "\"" or $c eq "\'" ) { 
      if( $c eq $close ) { 
        $close = ""; 
      } elsif( not $close ) {
        $close = "$c";
      } else {
        $new .= "$c";
      }
    } elsif( not $close ) {
      $new .= "$c";
    }
  }
  die "*** cannot parse string: $token\n" if $close;

  return $new;
}

sub elim_pars {

  my $token = shift;

  my $i = index($token,'(');
  return $token if $i == -1;

  my $new = substr($token,0,$i);
  my $rest = substr($token,$i);

  my $m = 0;
  my @f = split(//,$rest);
  foreach my $c (@f) {
    $m++ if $c eq '(';
    $new .= "$c" if $m == 0;
    $m-- if $c eq ')';
    die "*** cannot parse parenthesis: $token\n" if $m < 0;
  }
  die "*** cannot parse parenthesis: $token\n" if $m > 0;

  return $new;
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub handle_call {

  my ($self,$from,$name) = @_;

  $name = $self->upperlow($name);

  #print STDERR "    ........... call from $from to $name ($type)\n";

  $self->{called}->{$name}++;
  unless ( $self->{calling}->{$from} and 
		$self->{calling}->{$from} =~ / $name / ) {
    $self->{calling}->{$from} .= " $name ";
  }
  unless ( $self->{calledby}->{$name} and 
		$self->{calledby}->{$name} =~ / $from / ) {
    $self->{calledby}->{$name} .= " $from ";
  }
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub set_ignore_files {

	my ($self,$ignore) = @_;

	my %ignore_files = ();

	foreach my $file (@$ignore) {
          my $name = $self->upperlow($file);
	  $ignore_files{ $name } = 1;
	}

	$self->{ignore_files} = \%ignore_files;
}

sub set_ignore_routines {

	my ($self,$ignore) = @_;

	my %ignore_routines = ();

	foreach my $file (@$ignore) {
          my $name = $self->upperlow($file);
	  $ignore_routines{ $name } = 1;
	}

	$self->{ignore_routines} = \%ignore_routines;
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub calling_sequence {

  my ($self,$name,$ignore) = @_;

  return if( $ignore and $ignore->{$name} );		#programs to ignore

  print STDERR "calling sequence for $name:\n" unless $self->{level};

  $self->push_calling_stack($name);
  $self->{calling_times}->{$name}++;

  my $def = $self->{defined}->{$name};
  my $is_defined = "";
  $is_defined = "(-)" unless $def;

  if( $self->{maxlevel} and $self->{level} > $self->{maxlevel} ) {
    print STDERR "******* maxlevel: $name\n";
    $self->{maxlevel_reached}++;
    #$self->print_calling_stack(); exit 1;
    return;
  }
  $self->{level}++;

  $self->{maxl} = $self->{level} if $self->{maxl} < $self->{level};

    my $oldspace = $self->{space};
    $self->{space} .= "   ";
    print "$oldspace $name $is_defined\n";

    my $list = trim($self->{calling}->{$name});
    if( $list ) {
      my @list = split("  ",$list);
      foreach my $sub (@list) {
        $self->calling_sequence($sub);
      }
    }
    $self->{space} = $oldspace;

  $self->pop_calling_stack($name);

  $self->{level}--;
}

#----------------------------------------------------------------

sub print_calling_times {

  my ($self,$mincount) = @_;

  return if $mincount < 1;

  my $ct = $self->{calling_times};

  my @sorted = sort { $ct->{$a} <=> $ct->{$b} } keys %$ct;

  print "calling times (over $mincount):\n";

  foreach my $name (@sorted) {
    my $count = $ct->{$name};
    next if $count < $mincount;
    print "$count   $name\n";
  }
}

#----------------------------------------------------------------

sub push_calling_stack {

  my ($self,$name) = @_;

  if( $self->{calling_names}->{$name} ) {
    print_calling_stack();
    die "recursive calling stack: $name\n";
  } else {
    $self->{calling_names}->{$name} = 1;
    push(@{$self->{calling_stack}},$name);
  }
}

sub pop_calling_stack {

  my ($self,$name) = @_;

  if( $self->{calling_names}->{$name} ) {
    $self->{calling_names}->{$name} = 0;
    pop(@{$self->{calling_stack}});
  } else {
    print_calling_stack();
    die "corrupt calling stack: $name\n";
  }
}

sub print_calling_stack {

  my ($self) = @_;

  my $line = join(" ",@{$self->{calling_stack}});
  print STDERR "$line\n";
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub set_changed {

  my ($self,$file) = @_;

  my $fitem = $self->{files}->{$file};
  $fitem->{changed} = 1;
}

sub write_files {

  my ($self,$quiet,$all) = @_;

  foreach my $filename (keys %{$self->{files}}) {
    #print STDERR "calling write_file for $filename\n";
    $self->write_file($filename,$quiet,$all)
  }
}

sub write_file {

  my ($self,$filename,$quiet,$all) = @_;

  my $new_filename = "$filename.new";

  my $fitem = $self->{files}->{$filename};
  my $sequence = $fitem->{sequence};
  my $changed = $fitem->{changed};

  return unless $all or $changed;
  
  print STDERR "writing file $filename -> $new_filename\n" unless $quiet;
  open(NEW,">$new_filename");

  foreach my $ritem (@$sequence) {
    my $name = $ritem->{name};
    my $code = $ritem->{code};
    foreach my $line (@$code) {
      print NEW "$line\n";
    }
  }

  close(NEW);
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub upperlow {
  my ($self,$name) = @_;
  if( $self->{upper} ) {
    return uc($name);
  } else {
    return lc($name);
  }
}

sub trim {

  my $line = shift;

  return "" unless $line;

  $line =~ s/^\s+//;
  $line =~ s/\s+$//;

  return $line;
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
1;
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
