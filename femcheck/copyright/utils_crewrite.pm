#!/usr/bin/perl -ws

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------


#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

sub handle_c_revlog
{
  my ($rev,$header0,$header1,$header2) = @_;

  make_month_dates();

  if( $::cstyle_revlog ) {
    $header0 = clean_c_header($header0);
    $header1 = clean_c_header($header1);
    $header2 = clean_c_header($header2);
    $rev = rewrite_c_revlog($rev)
  } elsif( @$rev ) {
    print "*** file $::file already has new revison log... not rewriting\n";
  }

  return ($rev,$header0,$header1,$header2);
}

sub clean_c_header
{
  my $text = shift;

  foreach (@$text) {
    next unless $_;
    s/\s*\*\s*$//;
    $_ = " *" unless $_;	#deleted initial star... put it back
  }

  return $text;
}

sub rewrite_c_revlog
{
  my $rev = shift;

  my $nrev = @$rev;
  return $rev unless $nrev;
  if( $nrev == 1 ) {
    print STDERR "incomplete revision log in file $::file\n";
    return $rev;
  }
  my @aux = @$rev;
  my $first = shift @aux;
  unless( $first =~ /Revision History:/ ) {
    print STDERR "cannot parse revision header in file $::file\n";
    print STDERR "$first\n";
    return $rev;
  }
  my $second = shift @aux;
  unless( $second =~ /^\s*\*\s+\*\s*$/ ) {
    unshift(@aux,$second);
  }

  my $error = 0;
  @aux = reverse(@aux);
  my @new = ();
  my @conti = ();
  push(@new," * revision log :");
  push(@new," *");
  foreach (@aux) {
    #print "parsing $_\n";
    s/\s*\*\s*$//;
    #if( /^\s*\*\s+(\S)\s+(.+)\*\s*$/ ) {
    if( /^\s*\*\s+(\.\.\.)\s+(.+)$/ ) {
      my $text = $2;
      my $line = " * ...\t\tggu\t$text";
      unshift(@conti,$line);
    } elsif( /^\s*\*\s+(\S+)\s+(.+)$/ ) {
      my $date = $1;
      my $text = $2;
      $date = translate_c_date($date);
      my $line = " * $date\tggu\t$text";
      push(@new,$line);
      if( @conti ) {
        push(@new,@conti);
	@conti = ();
      }
    } else {
      print STDERR "cannot parse revision log in file $::file\n";
      print STDERR "$_\n";
      $error++;
    }
  }

  if( $error ) {
    return $rev;
  } else {
    return \@new;
  }
}

sub translate_c_date
{
  my $date = shift;

  my ($day,$month,$year);

  if( $date =~ /^(\d\d)-(\S+)-(\d\d\d\d):*$/ ) {
    $day = $1;
    $month = $2;
    $year = $3;
  } elsif( $date =~ /^(\S\S)-(\S+)-(\d{2,4}):*$/ ) {
    $day = $1;
    $month = $2;
    $year = $3;
    $year += 2000 if $year < 50;
    $year += 1900 if $year < 1900;
    $day = "01" if $day eq "..";
    if( $day < 1 or $day > 31 ) {
      print "*** error parsing day: $date ($::file)\n";
      return $date;
    }
    $month = "Jan" if $month eq "...";
  } else {
    print "*** error parsing date: $date ($::file)\n";
    return $date;
  }

  $month = $::month_dates{$month};
  return $date unless $month;

  unless( $month ) {
    print "*** error parsing month: $date ($::file)\n";
    return $date;
  }

  $date="$day.$month.$year";
  return $date;
}

sub revadjust_for_c
{
  print STDERR "   adjusting revlog for c\n";

  foreach (@_) {
    s/^./ */;
  }
  pop(@_);
  unshift(@_," *");
  #unshift(@_,"");
  unshift(@_,$::divisor_start);
  push(@_," *");
  push(@_,$::divisor_end);
  push(@_,"");

  return @_
}

sub make_month_dates {

  %::month_dates = (
         'Jan' => '01'
        ,'Feb' => '02'
        ,'Mar' => '03'
        ,'Apr' => '04'
        ,'May' => '05'
        ,'Jun' => '06'
        ,'Jul' => '07'
        ,'Aug' => '08'
        ,'Sep' => '09'
        ,'Oct' => '10'
        ,'Nov' => '11'
        ,'Dec' => '12'
        ,''    => ''
    );
}

#-------------------------------
1;
#-------------------------------


