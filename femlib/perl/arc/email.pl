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
# pleas use now address.pm

#$email_pattern = '[\w-%.]+@[\w-.]+\.[a-zA-Z]{2,3}';
$email_pattern = '[\w-%.]+@[\w-.]+';
$email_pattern1 = '[\w-.]+';			#just name

#---------------------------------------------------------------------

sub get_address {
  my $address = shift;
  my ($email,$name) = parse_address($address);
  return ($email,$name);
}

sub get_name {
  my $address = shift;
  my ($email,$name) = parse_address($address);
  return $name;
}

sub get_email {
  my $address = shift;
  my ($email,$name) = parse_address($address);
  return $email;
}

sub get_email_name {
  my $address = shift;
  my ($email,$name) = parse_address($address);
  $email =~ s/^(.*)\@(.*)$/$1/;
  return $email;
}

sub get_email_host {
  my $address = shift;
  my ($email,$name) = parse_address($address);
  $email =~ s/^(.*)\@(.*)$/$2/;
  return $email;
}

#---------------------------------------------------------------------

sub parse_address {

  my $address = shift;

  #chomp($address);

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

  $email = uniform_email( $email );
  $name = uniform_name( $name );

  return ($email,$name);
}

#---------------------------------------------------------------------

sub uniform_email {

  my $email = shift;

  return lc($email);
}

sub uniform_name {

  my $name = shift;

  $name =~ s/\"//g;
  $name =~ s/[\"()]//g;
  $name =~ s/,.*//g;
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

