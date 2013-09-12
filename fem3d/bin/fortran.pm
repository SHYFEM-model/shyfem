#!/usr/bin/perl -ws
#
# parses fortran routines
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
# $fortran->check_uniquness();
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
				 defined	  =>	{}
				,functions	=>	{}
				,called		  =>	{}
				,calling	  =>	{}
				,calledby	  =>	{}
				,routines	  =>	{}

				,ignore_files	      =>	{}
				,ignore_routines	  =>	{}

				,upper	    =>	0

				,file		    =>	""
				,line		    =>	""
				,insection	=>	""
				,code		    =>	[]

				,calling_names		    =>	{}
				,calling_stack		    =>	[]
			};

	bless $self;
	return $self;
}

#%defined = ();		# how often is this routine defined
#%functions = ();		# how often is this function defined
#%called = ();		# how often is this routine called
#%calling = ();		# routine calls these routines
#%calledby = ();		# routine is called by these routines
#%routines = ();		# complete info on routines

#$upper             # convert to upper case, else use lower case
#$code = "";
#$insection = "";
#$line = "";
#$file = "";
#$space = "";
#$level = 0;
#$maxl = 0;			# maximum level reached
#$maxlevel_reached = 0;

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub parse_files {

 my ($self,$quiet) = @_;

 $quiet = 0 unless defined $quiet;

 my $ignore_files = $self->{ignore_files};
 my $ignore = 0;

 while (<>) {

  if( $self->{file} ne $ARGV ) {
    print STDERR "reading $ARGV\n" unless $quiet;
    $self->{file} = $ARGV;
    $self->{line} = 0;
		my $file_name = $self->upperlow($ARGV);
		$ignore = $self->{ignore_files}->{$file_name};
		print STDERR "ignoring file $ARGV\n" if $ignore;
  }

	next if $ignore;

  chomp;

  $self->{line}++;

  my $orig = $_;
  s/^\s*\d+\s+/ /;	# get rid of label number
  s/\s+//g;		# no space
  s/\!.*//g;		# line comment

  if( /^END$/i ) { $self->insert_item("E",""); }
  if( /^END(FUNCTION|SUBROUTINE)(\w*)/i ) { $self->insert_item("E",$2); }
  if( /^PROGRAM(\w+)/i ) { $self->insert_item("P",$1); }
  if( /^SUBROUTINE(\w+)/i ) { $self->insert_item("S",$1); }
  if( /^FUNCTION(\w+)/i ) { $self->insert_item("F",$1); }
  if( /^BLOCKDATA(\w+)/i ) { $self->insert_item("B",$1); }
  if( /^(REAL|INTEGER|DOUBLEPRECISION|LOGICAL)FUNCTION(\w+)/i ) { 
    $self->insert_item("F",$2); 
  }

  push(@{$self->{code}},$orig) if $self->{code};
 }
}

sub insert_item {

  my ($self,$type,$name) = @_;

  $name = $self->upperlow($name);

  if( $type eq 'E' ) {		# close section
    unless($self->{insection}) {
      print STDERR "file: $self->{file}   line: $self->{line}\n";
      die "end statement found but in no section...\n";
    }
    $name = $self->{insection};
    my $item = $self->{routines}->{$name};
    $item->{code} = $self->{code};
    $self->{code} = [];
    $self->{insection} = "";
  } else {			# open new section
    if($self->{insection}) {
      print STDERR "file: $self->{file}   line: $self->{line}\n";
      die "routine found but already in section: $self->{insection}\n";
    }
    $self->{code} = [];
    my $item = {};
    $item->{name} = $name;
    $item->{type} = $type;
    $item->{file} = $self->{file};

    $self->{routines}->{$name} = $item;
    $self->{defined}->{$name}++;
    $self->{functions}->{$name}++ if $type eq "F";
    $self->{insection} = $name;
  }
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub show_functions {

  my ($self) = @_;

  print STDERR "list of functions:\n";
  foreach my $name (keys %{$self->{functions}}) {
    my $count = $self->{defined}->{$name};
    print STDERR "$count   $name\n";
  }
}

sub check_uniquness {

  my ($self) = @_;

  foreach my $name (keys %{$self->{defined}}) {
    my $count = $self->{defined}->{$name};
    if( $count != 1 ) {
      print STDERR "multiple definition of routine: $count  $name\n";
    }
  }
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub parse_routines {

  my ($self) = @_;

  foreach my $name (keys %{$self->{defined}}) {
    my $item = $self->{routines}->{$name};
    $self->parse_routine($item);
  }
}

sub parse_routine {

  my ($self,$item) = @_;

  my $name = $item->{name};
  my $type = $item->{type};
  my $code = $item->{code};

  #print STDERR "parsing $name ($type)\n";

  unless( $code ) {
    print STDERR "no code for routine $name\n";
    return;
  }

  $self->{insection} = $name;
  my $nline = 0;

  foreach my $line (@$code) {
    $nline++;
    next if $nline == 1;	# skip definition
    $line =~ s/^\s*\d+\s+/ /;	# get rid of label number

    next if $line =~ /^[cC!]/;		# get rid of line comment

    if( $line =~ /^\s+CALL\s+(\w+)/i ) { $self->handle_call("C",$1); }
    if( $line =~ /^\s+IF\s*\(.+\)\s*CALL\s+(\w+)/i ) { 
	$self->handle_call("C",$1); 
    }
    while( $line =~ /(\w+)\s*\(/ig ) { 
      my $func = $self->upperlow($1);
      if( $self->{functions}->{$func} ) { 
	#print STDERR "function call: $func in routine $name\n";
        $self->handle_call("F",$func);
      }
    }

  }
}

sub handle_call {

  my ($self,$type,$name) = @_;

  $name = $self->upperlow($name);
  my $insection = $self->{insection};

  #print STDERR "    ........... call from $insection to $name ($type)\n";

  $self->{called}->{$name}++;
  unless ( $self->{calling}->{$insection} and 
		$self->{calling}->{$insection} =~ / $name / ) {
    $self->{calling}->{$insection} .= " $name ";
  }
  unless ( $self->{calledby}->{$name} and 
		$self->{calledby}->{$name} =~ / $insection / ) {
    $self->{calledby}->{$name} .= " $insection ";
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

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub calling_sequence {

  my ($self,$name,$ignore) = @_;

  return if( $ignore and $ignore->{$name} );		#programs to ignore

  $self->push_calling_stack($name);
  $::calling_times{$name}++;

  my $def = $::defined{$name};
  my $is_defined = "";
  $is_defined = "(-)" unless $def;

  if( $::maxlevel and $::level > $::maxlevel ) {
    print STDERR "******* maxlevel: $name\n";
    $::maxlevel_reached++;
    #$self->print_calling_stack(); exit 1;
    return;
  }
  $::level++;

  $::maxl = $::level if $::maxl < $::level;

    print "$::space $name $is_defined\n";
    my $oldspace = $::space;
    $::space .= "   ";
    my $list = $::calling{$name};
    if( $list ) {
      $list =~ s/,$//;
      my @list = split(",",$list);
      foreach my $sub (@list) {
        $self->calling_sequence($sub);
      }
    }
    $::space = $oldspace;

  $self->pop_calling_stack($name);

  $::level--;
}

#----------------------------------------------------------------

sub print_calling_times {

  return if $::mincount < 1;

  my @sorted = sort by_count keys %::calling_times;
  #my @sorted = sort keys %::calling_times;

  my $i = 0;
  foreach my $key (@sorted) {
    my $count = $::calling_times{$key};
    next if $count < $::mincount;
    print "calling times (over $::mincount):\n" if $i == 0;
    print "$key   $count\n";
    $i++;
  }
}

#----------------------------------------------------------------
#----------------------------------------------------------------
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

sub upperlow {
  my ($self,$name) = @_;
  if( $self->{upper} ) {
    return uc($name);
  } else {
    return lc($name);
  }
}

sub by_count {

  if( $::calling_times{$a} < $::calling_times{$b} ) {
    return -1;
  } elsif( $::calling_times{$a} > $::calling_times{$b} ) {
    return +1;
  } else {
    return 0;
  }
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
1;
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
