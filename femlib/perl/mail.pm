#!/usr/bin/perl -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# mail utility routines
#
# Usage:
#
# #!/usr/bin/env perl
#
# use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
# use mail;
# 
# $file = $ARGV[0];
# 
# my $mbox = new mail;
# $mbox->open($file);
# $mbox->info;
#
# or more compact
#
# my $mbox = new mail($file);
# $mbox->info;
#
#
##############################################################

use strict;
use message;

package mail;

##############################################################

sub new
{
    my ($pck,$folder) = @_;

    my $self;

    $self =	{
	    		 folder		=>	undef
			,fhandle	=>	undef
			,msg		=>	[]
			,nmsg		=>	0
			,imsg		=>	0
			,is_mbox	=>	-1
			,envline	=>	undef
			,headline	=>	undef
		};

    bless $self;

    if( $folder ) {
      $self->open($folder);
    }

    return $self;
}

##############################################################

sub open
{
    my ($self,$folder) = @_;

    #print STDERR "opening folder $folder\n";
    my $imsg = 0;
    my $msgs = $self->{msg};

    if( $folder ne "-" ) {
      open(FOLDER,"<$folder") or die "Cannot open folder: $folder\n";
      $self->{fhandle} = \*FOLDER;
    } else {
      $self->{fhandle} = \*STDIN;
      #print STDERR "reading from stdin...\n";
    }

    $self->{is_mbox} = $self->is_mbox();

    while( my $msg = $self->read_next_message() ) {
      push(@$msgs,$msg);
      $imsg++;
      #print STDERR "message read: $imsg\n";
    }
    $self->{nmsg} = $imsg;
    $self->{folder} = $folder;

    if( $folder ne "-" ) {
      close(FOLDER);
      $self->{fhandle} = undef;
    }
}

sub info	#retuns number of messages, -1 if not a mailbox
{
    my ($self) = @_;

    my $nmsg = 0;

    #print STDERR "folder: $self->{folder} with $self->{nmsg} messages";
    #print STDERR " ($self->{is_mbox})\n";

    if( $self->{is_mbox} ) {
      #print STDERR "folder: $self->{folder} with $self->{nmsg} messages\n";
      print STDERR "$self->{nmsg}    $self->{folder}\n";
      $nmsg = $self->{nmsg};
    } else {
      #print STDERR "folder: $self->{folder} is not a mailbox\n";
      print STDERR "not a mailbox    $self->{folder}\n";
      $nmsg = -1;
    }

    return $nmsg;
}

sub reset_messages
{
    my ($self) = @_;

    $self->{imsg} = 0;
}

sub next_message
{
    my ($self) = @_;

    my $msgs = $self->{msg};
    my $nmsg = $self->{nmsg};
    my $imsg = $self->{imsg};

    return undef if $imsg >= $nmsg;

    $self->{imsg} += 1;

    return $$msgs[$imsg];
}

sub get_message
{
    my ($self,$imsg) = @_;

    my $msgs = $self->{msg};
    my $nmsg = $self->{nmsg};

    return undef if $nmsg == 0;

    $imsg = $imsg % $nmsg;
    $self->{imsg} = $imsg;

    return $$msgs[$imsg];
}

##############################################################

sub read_next_message
{
    my ($self) = @_;

    return "" unless( $self->{is_mbox} );

    my $env = $self->read_envelop();
    return undef unless $env;

    my $head = $self->read_header();
    my $body = $self->read_body();

    my $msg = new message($env,$head,$body);

    return $msg;
}

sub read_envelop
{
    my ($self) = @_;

    my $envelop;

    if( $self->{envline} ) {
      $envelop = $self->{envline};
      $self->{envline} = undef;
    } else {
      my $FH = $self->{fhandle};
      $envelop = <$FH>;
    }

    return $envelop;
}

sub read_header
{
    my ($self) = @_;

    my @header = ();
    my $FH = $self->{fhandle};

    if( $self->{headline} ) {
      push(@header,$self->{headline});
      $self->{headline} = undef;
    }

    while( my $line = <$FH> ) {
      last if $line =~ /^\s*$/ ;	# look for empty line
      push(@header,$line);
    }

    return \@header;
}

sub read_body
{
    my ($self) = @_;

    my @body = ();
    my $FH = $self->{fhandle};

    while( my $line = <$FH> ) {
      if( $self->is_envelop($line) ) {
	my $first = <$FH>;
	if( $self->is_header_line($first) ) {
          $self->remember_envelop($line,$first);
          last;
	} else {	# false envelop
	  #print STDERR "*** fake envelop found ...\n";
	  #print STDERR "$line";
	  #print STDERR "$first";
          push(@body,$line);
          push(@body,$first);
	  next;
        }
      }
      push(@body,$line);
    }

    return \@body;
}

sub is_header_line
{
    my ($self,$line) = @_;

    #print STDERR "check header: $line\n";

    if( $line =~ /^\S+: / ) {
      return 1;
    } else {
      return 0;
    }
}

sub is_envelop
{
    my ($self,$line) = @_;

    if( $line and $line =~ /^From / ) {
      return 1;
    } else {
      return 0;
    }
}

sub remember_envelop
{
    my ($self,$line,$first) = @_;

    $self->{envline} = $line;
    $self->{headline} = $first;
}

sub is_mbox
{
    my ($self) = @_;

    my $FH = $self->{fhandle};

    if( $self->{is_mbox} == -1 ) {
      my $line = <$FH>;
      #if( not defined $line ) {
      #  $self->{is_mbox} = 1;		#empty
      #  return $self->{is_mbox};
      #}
      $self->{is_mbox} = $self->is_envelop($line);
      $self->remember_envelop($line);
     
    }

    return $self->{is_mbox};
}

##############################################################


###################################
#&test_str(@ARGV);
###################################
1;
###################################

