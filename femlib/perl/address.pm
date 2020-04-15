#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# example of usage:
#
# use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
#
# use address;
#
# my $address = new address;
# $address->set_address($from);
# $name = $address->get_name();
#
# or
#
# my $address = new address($from);
# $name = $address->get_name();
#
#-----------------
#
# example:
#
# $from = "Georg U.  <georg@aol.com>";
# my $address = new address($from);
#
# my ($email,$name) = $address->get_address()	# ("georg@aol.com","Georg U.")
# my $name = $address->get_name();		# "Georg U."
# my $email = $address->get_email();		# "georg@aol.com"
# my $ename = $address->get_email_name();	# "georg"
# my $ehost = $address->get_email_host();	# "aol.com"
#
# for multiple addresses the above routines only return the first address
# to get all addresses use
#
# my ($count,$email,$name) = $address->get_addresses()
#
# where $email,$name are array references and $count is the number of 
# addresses found
#
#-----------------

use strict;

package address;

#---------------------------------------------------------------------

sub new {

    my ($self,$address) = @_;

    $self =     {
                         address           =>      undef
                };

    bless $self;
    $self->set_address($address) if $address;
    return $self;
}

#---------------------------------------------------------------------

sub set_address {
  my ($self,$address) = @_;
  $self->{address} = $address;
}

sub get_address {
  my $self = shift;
  my ($email,$name) = $self->parse_address();
  return ($email,$name);
}

sub get_name {
  my $self = shift;
  my ($email,$name) = $self->parse_address();
  return $name;
}

sub get_email {
  my $self = shift;
  my ($email,$name) = $self->parse_address();
  return $email;
}

sub get_email_name {
  my $self = shift;
  my ($email,$name) = $self->parse_address();
  $email =~ s/^(.*)\@(.*)$/$1/;
  return $email;
}

sub get_email_host {
  my $self = shift;
  my ($email,$name) = $self->parse_address();
  $email =~ s/^(.*)\@(.*)$/$2/;
  return $email;
}

#---------------------------------------------------------------------

sub get_addresses {

  my $self = shift;

  return $self->parse_addresses();
}

sub parse_addresses {

  my $self = shift;

  my $i = 0;
  my @email = ();
  my @name = ();

  my $addresses = $self->parse_multiple_addresses();

  foreach my $address (@$addresses) {
    my ($email,$name) = $self->parse_single_address($address);
    push(@email,$email);
    push(@name,$name);
    $i++;
  }

  return ($i,\@email,\@name);
}

sub parse_address {

  my $self = shift;

  my $addresses = $self->parse_multiple_addresses();
  my $address = $addresses->[0];

  #$address = $self->{address};	#old way
  return $self->parse_single_address($address);
}

sub parse_single_address {

  my ($self,$address) = @_;

  my $email_pattern = '[\w-%.]+@[\w-.]+';
  my $email_pattern1 = '[\w-.]+';			#just name

  my ($email,$name);

  $address =~ s/^\s+//;
  $address =~ s/\s+$//;
  $address =~ s/\s+/ /g;

  return unless $address;

  if( $address =~ /<($email_pattern)>/ ) {                  # name <email>
    $email = $1;
    $name = $address;
    $name =~ s/<$email>//;
  } elsif( $address =~ /<($email_pattern1)>/ ) {            # name <email>
    $email = $1;
    $name = $address;
    $name =~ s/<$email>//;
  } elsif( $address =~ /\[SMTP:($email_pattern)\]/ ) {      # name [SMTP:email]
    $email = $1;
    $name = $address;
    $name =~ s/\[SMTP:$email\]//;
  } elsif( $address =~ /^($email_pattern)\s+\((.*)\)$/ ) {  # email (name)
    $email = $1;
    $name = $2;
  } elsif( $address =~ /^($email_pattern1)\s+\((.*)\)$/ ) {  # email (name)
    $email = $1;
    $name = $2;
  } elsif( $address =~ /^($email_pattern)$/ ) {             # email
    $email = $1;
    $name = "";
  } elsif( $address =~ /^($email_pattern1)$/ ) {            # email
    $email = $1;
    $name = "";
  } elsif( $address =~ /Recipient list suppressed/i ) {
    $email = "";
    $name = "";
  } elsif( $address =~ /recipient list not shown/i ) {
    $email = "";
    $name = "";
  } elsif( $address =~ /Undisclosed[ -]Recipient/i ) {
    $email = "";
    $name = "";
  } elsif( $address =~ /MAILER-DAEMON/i ) {
    $email = "";
    $name = "";
  } else {
    print STDERR "*** Cannot parse address: $address\n";
  }

  #print "parse: |$address|$email|$name|\n";

  $email = $self->uniform_email( $email );
  $name  = $self->uniform_name ( $name  );

  return ($email,$name);
}

sub parse_multiple_addresses {

  my $self = shift;
  my $address = $self->{address};

  my $orig = $address;
  my @addresses = ();
  my ($act,$div);
  my $debug = 0;

  $address =~ s/\n/ /g;

  while( $address ) {
	#print "|$address|\n";
    if( $address =~ /(.*?)([\"\(,])(.*)/ ) {
      $act .= $1;
      $div = $2;
      $address = $3;
      if( $div eq "," ) {
	push(@addresses,$act);
	$act = "";
      } else {
	if( $div eq '"' and $address =~ /(.*?)\"(.*)/ ) {
	  $act .= '"' . $1 . '"';
	  $address = $2;
	} elsif( $div eq '(' and $address =~ /(.*?)\)(.*)/ ) {
	  $act .= '(' . $1 . ')';
	  $address = $2;
	} else {
	  print STDERR "*** Cannot parse address: $orig ($address)\n";
	  return [];
	}
      }
    } else {
      $act .= $address;
      $address = "";
    }
  }

  push(@addresses,$act) if $act;
  my $n = @addresses;

  if( $debug and $n > 1 ) {
	print STDERR "** multiple addresses found: $n\n";
	print STDERR "      $orig\n";
	foreach my $address (@addresses) {
	  print STDERR "        $address\n";
	}
  }

  return \@addresses;
}

#---------------------------------------------------------------------

sub uniform_email {

  my ($self,$email) = @_;

  return lc($email);
}

sub uniform_name {

  my ($self,$name) = @_;

  $name =~ s/\"//g;
  $name =~ s/[\"()]//g;
  $name =~ s/^\s*\'(.*)\'\s*$/$1/g;	# 'georg umgiesser'
  #$name =~ s/,.*//g;
  $name =~ s/\s*Dr\.\s*//g;
  $name =~ s/\s+by way of.*$//ig;
  $name =~ s/\s+tramite.*$//ig;

  $name =~ s/^\s+//;
  $name =~ s/\s+$//;

  return $name;
}

#---------------------------------------------------------------------
1;
#---------------------------------------------------------------------

