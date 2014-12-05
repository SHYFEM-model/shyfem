#!/usr/bin/perl -ws
#
# parses fortran routines
#
# version 1.2		02.12.2014
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

				,routines	    =>	{}
				,files		    =>	{}

				,ignore_files	    =>	{}
				,ignore_routines    =>	{}

				,upper		    =>	0
				,no_main	    =>	0

				,comnum		    =>	0
				,implicit	    =>	0

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
    $name = $self->make_comment_name();
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

sub make_comment_name {

  my ($self) = @_;

  $self->{comnum}++;
  return "comment_" . $self->{comnum};
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
  $ritem->{code} = $code;

  my $fitem = $self->{files}->{$file};
  my $seq = $fitem->{sequence};
  push(@$seq,$ritem);

  return if $type eq "comment";
  $self->{functions}->{$name}++ if $type eq "function";

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

  #print STDERR "reading routine: $name $type $next\n";

  while (<FILE>) {

    chomp;
    $self->{nlines}++;
    $orig = $_;
    push(@code,$_);

    s/^\s*\d+\s+/ /;	# get rid of label number
    s/\s+//g;		# no space
    s/\!.*//g;		# line comment

    if( /^INTERFACE$/i ) { 
      $in_interface = 1;
    } elsif( /^ENDINTERFACE$/i ) { 
      $in_interface = 0;
    }
    next if $in_interface;
    
    last if( /^END$type$name/i );
    last if( /^END$type/i );
    last if( /^END$/i );
  }

  die "*** cannot find end statement: $orig\n" unless ( /^END/i );

  return ($name,$type,\@code);
}

sub parse_initial {

  my ($self,$next) = @_;

  $_ = $next;
  s/^\s*\d+\s+/ /;	# get rid of label number
  s/\s+//g;		# no space
  s/\!.*//g;		# line comment

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
    foreach my $sitem (@$sequence) {
      my $name = $sitem->{name};
      print "    $name\n";
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

  my $implicit = 0;

  #print STDERR "parsing $name ($type)\n";

  unless( $code ) {
    print STDERR "*** no code for routine $name\n";
    return;
  }

  my $nline = 0;

  foreach my $cline (@$code) {
    $nline++;
    next if $nline == 1;	# skip definition
    my $line = $cline;

    $line =~ s/^\s*\d+\s+/ /;	# get rid of label number
    next if $line =~ /^[cC!]/;		# get rid of line comment

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

    $implicit = 1 if $line =~ /^\s*implicit\s*none/;
  }

  if( $self->{implicit} and not $implicit ) {
    print STDERR "no implicit none in $name (file $file)\n";
    my $fitem = $self->{files}->{$file};
    $fitem->{implicit} .= " $name ";
  }
  $ritem->{implicit} = $implicit;
}

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

  foreach my $sitem (@$sequence) {
    my $name = $sitem->{name};
    my $code = $sitem->{code};
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
