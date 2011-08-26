#!/usr/bin/perl -s

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use mail;
use message;
use address;

use strict;

package spam;

#-------------------------------------------------------------

sub new {

    my $self;

    $self =     {
                         spam_list           =>      undef
		};

    $self->{spam_list} = read_spam_list();

    bless $self;
    return $self;
}

#-------------------------------------------------------------

sub is_spam {

  my ($self,$msg) = @_;

  my $spam_list = $self->{spam_list};

  my $from = $msg->get_header("From");
  my $address = new address($from);
  my $email = $address->get_email();
  my $host = $address->get_email_host();

  if( $$spam_list{$host} ) {
    return 1;				#spam
  } else {
    return 0;				#no spam
  }
}

#-------------------------------------------------------------

sub read_spam_list {

  my %spam = ();
  my $spamfile = "/home/georg/Mail/proctest/SPAM";

  open(FILE,"<$spamfile") or die "Cannot read file: $spamfile\n";

  while( my $line = <FILE>) {

    chomp($line);
    next if $line =~ /^\s*$/;		#nothing on line
    next if $line =~ /^\s*#/;		#comment

    #print STDERR "spam... $line\n";
    $spam{$line} = 1;

  }

  close(FILE);

  return \%spam;
}

#-------------------------------------------------------------
1;
#-------------------------------------------------------------

