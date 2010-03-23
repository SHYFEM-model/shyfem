#!/usr/bin/perl -w
#
# Mailbox routines
#
#########################################

use strict;
use GguFile;

#########################################

package ForFile;
 
sub new
{
    my $self;

    my @routines = ();

    $self = 	{
			 file		=>	undef
			,ggufile	=>	undef
			,lastline	=>	undef
			,act_routine	=>	undef
			,routines	=>	\@routines
		};

    bless $self;
    return $self;
}

sub open
{
    my ($self,$file) = @_;

    my $ggufile = new GguFile;
    $ggufile->open_file( $file );

    $self->{file} = $file;
    $self->{ggufile} = $ggufile;
}

sub close
{
    my $self = shift;

    my $ggufile = $self->{ggufile};

    $ggufile->close_file;
}

sub get_act_routine
{
    my $self = shift;

    return $self->{act_routine};
}

sub set_act_routine
{
    my ($self,$act_routine) = @_;

    if( $act_routine ) {
        $self->{act_routine} = $act_routine;
    } else {
        $self->{act_routine} = undef;
    }
}

sub close_act_routine
{
    my $self = shift;

    if( $self->{act_routine} ) {
	my $routines = $self->{routines};
        push(@$routines,$self->{act_routine});
    }
    $self->{act_routine} = undef;
}

sub get_next_statement
{
    my $self = shift;

    my $pline = $self->get_line;
    return "" unless $pline;

    if( &is_conti($pline) ) {
        die "impossible continuation line:\n$pline->{line}\n";
    }

    my $state = $pline;

    while( $pline = $self->get_line ) {
        if( &is_conti($pline) ) {
            $state->{line}      .=        $pline->{line};
            $state->{statement} .= " "  . $pline->{statement};
            $state->{comment}   .= "\n" . $pline->{comment};
            $state->{conti}++;
        } else {
            $self->unget_line($pline);
            last;
        }
    }

  $self->{act_statement} = $state;
  #print "new: $state->{line}\n";
  return $state;

}

sub get_line
{
    my $self = shift;

    my $ggufile = $self->{ggufile};
    my $pline;

    if( $self->{lastline} ) {
        $pline = $self->{lastline};
        $self->{lastline} = undef;
    } else {
        my $line = $ggufile->get_line;
        return "" unless $line;			#end of file
        $pline = &pre_parse($line);
        $pline->{line_number} = $ggufile->lines_read;
    }

    return $pline;
}

sub unget_line
{
    my ($self,$pline) = @_;

    $self->{lastline} = $pline;
}

sub pre_parse
{

  # this routine sets the following fields:
  # line, statement, conti, nstate, comment
  # in case of only comment it also sets keyword=comment
  # in case of empty line it also sets keyword=comment

  $_ = shift;

  return "" unless $_;

  my %pline = {};

  $pline{"line"} = $_;

  chomp;
  &untab;

  tr/A-Z/a-z/;		#make lowercase

  # treat comment line

  if( /^c(.*)/i || /^\*(.*)/ || /^\s*!(.*)/ ) {		#only comment
    $pline{comment} = $1;
    $pline{keyword} = "comment";
    return \%pline;
  } elsif( /^\s*$/ ) {					#empty line
    $pline{keyword} = "comment";
    $pline{comment} = "";
    return \%pline;
  }

  # treat trailing comment

  if( /^(.*)!(.*)/ ) {
	my $code = $1;
	my $rest = $2;
	unless( $rest =~ /\'/ ) {	#FIXME
		$_ = $code;
	}
        $pline{comment} = $rest;
  }

  # treat conti line and statement number

  my $n = substr($_,0,5);
  my $c = substr($_,5,1);
  my $rest = substr($_,6);

  if( $c && $c ne " " && $c ne "0" ) {
    $pline{conti} = 1;
  }

  $n =~ s/\s+//g;		#remove blanks
  if( $n =~ /(\d+)/ ) {
    $pline{nstate} = $1;
  }

  if( $pline{conti} && $pline{nstate} ) {
    die "Statement line on continuation number impossible:\n$_\n";
  }

  # treat rest

#FIXME next
  #$rest = &ggu_strings::subst_strings( $rest );		# "..." -> $STRING1$

  $rest =~ s/^\s*(.*)\s*$/$1/;		#remove blanks
  $rest =~ s/\s+//g;			#remove all blanks
  $pline{statement} = $rest;

  return \%pline;
}

#-------------------------------------------------------------------

sub is_conti {

  my $pline = shift;

  return $pline->{conti};
}

  
#-------------------------------------------------------------------

sub lower_case {

  my $string = shift;

  $string =~ tr/A-Z/a-z/;		#make lowercase

  return $string;
}

  
#-------------------------------------------------------------------

sub untab {

  my $tablength = 8;

  while( (my $i = index($_,"\t")) != -1 ) {
    my $j = $i % $tablength;
    $j = $tablength - $j;
    if( $j > 0 ) {
	substr($_,$i,1) = substr('           ',0,$j);
    }
  }
}
  
#-------------------------------------------------------------------

1;

#-------------------------------------------------------------------

