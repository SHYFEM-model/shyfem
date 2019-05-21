#!/usr/bin/perl -ws

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use strict;

sub get_revision_log
{

  my $in_revision = 0;
  my @revlog = ();

  while(<>) {

    chomp;

    check_copyright();

    if( $in_revision == 0 ) {
      if( /revision log :/) {		#start of revision log
	if( $::has_revision_log ) {
	  print "more than one revision log in file $ARGV\n";
	}
        $in_revision = 1;
        $::has_revision_log = 1;
      } else {
        check_revision(0);
      }
    } else {
      if( /^\s*$/ ) {			#empty line - end of revision log
        $in_revision = 0;
      } elsif( /^[!cC]\*\*\*/ ) {		#c*** line
        $in_revision = 0;
      } elsif( /^[!cC]\-\-\-/ ) {		#c--- line
        $in_revision = 0;
      } elsif( /^[!cC]\=\=\=/ ) {		#c=== line
        $in_revision = 0;
      } elsif( /^[!cC].*\s:\s*$/ ) {	#other comment line
        $in_revision = 0;
      } else {
        check_revision(1);
      }
      if( $in_revision ) {		#save revision log
        my $item = parse_revision_line($_);
        push(@revlog,$item) if $item;
      }
    }
  }
  return \@revlog;
}

sub is_start_of_revision_log
{
  my $line = shift;

  if( $line =~ /^(\S)\s+revision log :/ ) {	#start of revision log
    $::comment_char = $1;
    return 1;
  } else {
    return 0;
  }
}

sub is_end_of_revision_log
{
  my $line = shift;

  if( $line =~ /^\s*$/ ) {			#empty line
    return "$::final_revlog$line";
  } elsif( $line =~ /^[!cC]\*\*\*/ ) {		#c*** line
    return "$::final_revlog$line";
  } elsif( $line =~ /^[!cC]\-\-\-/ ) {		#c--- line
    return "$::final_revlog$line";
  } elsif( $line =~ /^[!cC]\=\=\=/ ) {		#c=== line
    return "$::final_revlog$line";
  } elsif( $line =~ /^[!cC].*\s:\s*$/ ) {	#other comment line
    return "$::final_revlog$line";
  } elsif( $line =~ /^[!cC]\s*$/ ) {		#empty comment line
    $::final_revlog .= $line;
    return 0;
  } else {
    $::final_revlog = "";
    return 0;
  }
}

sub parse_revision_line
{
    my $line = shift;

    return if $line =~ /^\S*$/;            #skip empty line
    return if $line =~ /^[cC!*]\s*$/;      #skip just comment line

    $line =~ s/^[cC!*]\s*//;            #strip comment in front

    my %item = ();
    if( $line =~ /^(\d\S+)\s+(\S+)\s+(.*)/ ) {
      $item{date} = $1;
      $item{name} = $2;
      $item{text} = $3;
      $item{idate} = make_idate($1);
    } elsif( $line =~ /^(\.\.\.)\s+(.*)/ ) {
      $item{date} = $1;
      $item{text} = $2;
      $item{idate} = -1;
    } else {
      die "cannot parse revlog (ARGV): $_\n";
    }

  return \%item;
}

sub read_revision_log
{
  my $file = shift;

  my @revlog = ();

  open(REV,"<$file") || die "Cannot open file: $file\n";
  while(<REV>) {
    chomp;
    next if /^[cC!]\s*$/;
    next if /^\s*$/;
    my $item = parse_revision_line($_);
    push(@revlog,$item) if $item;
  }
  close(REV);

  return \@revlog;
}

sub write_revision_log
{
  my ($file,$items) = @_;

  my $c = $::comment_char;
  $c = "!" unless $c;

  open(REV,">$file") || die "Cannot open file: $file\n";
  foreach my $item (@$items) {
    my $idate = $item->{idate};
    my $date = $item->{date};
    my $name = $item->{name};
    my $text = $item->{text};
    if( $idate == -1 ) {
      print REV "$c ...\t\t\t\t$text\n";
    } else {
      print REV "$c $date\t$name\t$text\n";
    }
  }
  close(REV);
}

sub check_dates_of_revision_log
{
  my ($ra,$file) = @_;

  $file = "" unless $file;

  my $idate_old = -1;

  foreach my $item (@$ra) {
    my $idate = $item->{idate};
    next if $idate == -1;
    if( $idate < $idate_old ) {
      print STDERR "error in revision log dates ($file): $idate $idate_old\n";
    }
    $idate_old = $idate;
  }
}

sub read_lines
{
  my $file = shift;

  my @lines = ();

  open(FILE,"<$file") || die "Cannot open file: $file\n";
  while(<FILE>) {
    push(@lines,$_);
  }
  close(FILE);

  return \@lines;
}

sub print_lines
{
  my $lines = shift;

  foreach my $line (@$lines) {
    print $line;
  }
}

sub show_lines
{
  my ($lines,$n,$text) = @_;

  my $l = @$lines;
  $n = 0 unless $n;

  print STDERR "show_line ($l): $text\n";

  if( $n == 0 or 2*$n > $l ) {
    foreach my $line (@$lines) {
      print STDERR $line;
    }
  } else {
    for( my $i=0; $i<$n; $i++ ) {
      print STDERR $lines->[$i];
    }
    print STDERR "...\n";
    for( my $i=$l-$n; $i<$l; $i++ ) {
      print STDERR $lines->[$i];
    }
  }
}

#--------------------------------------------------------------

sub check_copyright
{
  $::has_copyright = 1 if /^[cC!]\s+Copyright \(C\)/;
  if( /^[cC!]\s+This file is part of SHYFEM/ ) {
    $::has_shyfem = 1;
    $::is_manual = 1 if /\(m\)\s*$/;
  }
}

sub skip_over_copyright
{
  my $file = shift;

  my @copy = ();
  my @rest = ();

  open(FILE,"<$file") || die "Cannot open file: $file\n";

  while(<FILE>) {
    push(@copy,$_);
    check_copyright();
    if( /^!--------------------------/ ) {
      last if $::has_copyright;
    }
  }

  if( $_ ) {
    while(<FILE>) {
      push(@rest,$_);
    }
  }

  close(FILE);

  if( not $::has_copyright ) {
    print "file has no copyright...\n";
    return(\@copy,\@rest);
  } elsif( $::is_manual ) {
    print "file has manual copyright...\n";
    return(\@copy,\@rest);
  } else {
    return(\@copy,\@rest);
  }
}

sub skip_over_revision_log
{
  my $file = shift;

  my @copy = ();
  my @rest = ();
  my $in_revision = 0;

  open(FILE,"<$file") || die "Cannot open file: $file\n";

  while(<FILE>) {
    if( $in_revision ) {
      if( my $line = is_end_of_revision_log($_) ) {
        push(@rest,$line);
        last;
      }
      next;
    } else {
      if( is_start_of_revision_log($_) ) {
        $in_revision = 1;
        $::has_revision_log = 1;
        next;
      }
    }
    push(@copy,$_);
    check_copyright();
  }

  if( $_ ) {
    while(<FILE>) {
      push(@rest,$_);
    }
  }

  close(FILE);

  if( not $::has_copyright ) {
    print "file has no copyright...\n";
    return(\@copy,\@rest);
  } elsif( $::is_manual ) {
    print "file has manual copyright...\n";
    return(\@copy,\@rest);
  } else {
    return(\@copy,\@rest);
  }
}

#--------------------------------------------------------------

sub combine_revision_logs
{
  my ($file1,$file2) = @_;

  my $r1 = read_revision_log($file1);
  my $r2 = read_revision_log($file2);

  my @new = ();

  my $item1 = shift(@$r1);
  my $item2 = shift(@$r2);

  while(1) {
    if( not $item1 and not $item2 ) {
      last
    } elsif( not $item1 ) {
      push(@new,$item2);
      $item2 = shift(@$r2);
    } elsif( not $item2 ) {
      push(@new,$item1);
      $item1 = shift(@$r1);
    } elsif( $item1->{idate} <= $item2->{idate} ) {
      push(@new,$item1);
      if( $item2->{idate}-$item1->{idate} < 7 ) {	# skip git item
        $item2 = shift(@$r2);
      }
      $item1 = shift(@$r1);
    } else {
      push(@new,$item2);
      $item2 = shift(@$r2);
    }
  }

  return \@new;
}

sub combine_revision_log
{
  my ($file,$file2) = @_;

  my $r1 = read_revision_log($file);
  my $r2 = read_revision_log($file2);

  my %hash = ();

  foreach my $item (@$r2) {
    my $idate = $item->{idate};
    $hash{ $idate } = $item;
  }
  foreach my $item (@$r1) {
    my $idate = $item->{idate};
    $hash{ $idate } = $item;
  }

  my @keys = sort( keys( %hash ) );
  my @new = ();

  my $idate_old = -1;
  foreach my $key (@keys) {
    my $item = $hash{$key};
    my $idate = $item->{idate};
    my $diff = $idate - $idate_old;
    my $text = $item->{text};
    if( $diff < 10 and $text =~ /^changed/ ) {
      my $line = substr($text,0,22);
      print STDERR "  eliminated: $diff $idate_old $idate $line\n";
      next;
    }
    push(@new,$item);
    $idate_old = $idate;
  }

  return \@new;
}

#--------------------------------------------------------------

sub substitute_comment_char
{
  my $revlog = shift;

  my $c = $::comment_char;

  foreach (@$revlog) {
    s/^./$c/;
  }
}

sub substitute_revision_log
{
  my ($file,$revfile) = @_;

  my $revlog = read_lines($revfile);

  my ($copy,$rest) = skip_over_revision_log($file);
  my $c = $::comment_char;
  substitute_comment_char($revlog);

  print_lines($copy);

  if( $::has_copyright and not $::is_manual ) {
    print "$c revision log :\n$c\n";
    foreach (@$revlog) { print; }
    #print "\n";
  }

  print_lines($rest);
}

sub integrate_revision_log
{
  my ($file,$revfile) = @_;

  my $revlog = read_lines($revfile);

  my ($copy,$rest) = skip_over_copyright($file);

  print_lines($copy);

  if( $::has_copyright and not $::is_manual ) {
    print_revision_log($revlog);
  }

  print_lines($rest);
}

sub print_revision_log
{
  my $revlog = shift;

  print "\n";
  print "! revision log :\n";
  print "!\n";

  foreach (@$revlog) {
    print;
  }

  my $stars = '*'x74;

  print "\n";
  print "!$stars\n";
}

#--------------------------------------------------------------

sub check_revision {

  my $irv = shift;

  return if /^[!cC]\s*$/;
  return if /^[!cC]\*\*\*/;
  return if /^[!cC]\-\-\-/;

  if( my $iirv = check_new_revision() ) {
    if( $irv == 0 and $::warn ) {
      if( $iirv != 2 ) {	#not a continuation line
        print STDERR "new revision out of revision log ($ARGV): $_\n";
      }
    }
  } elsif( check_old_revision() ) {
    if( $irv == 0 and $::warn ) {
      print STDERR "old revision out of revision log ($ARGV): $_\n";
    }
  } elsif( $::obsolete and check_obsolete_revision() ) {
    ;
  } else {
    if( $irv ) {
      print STDERR "*** Cannot parse revision log ($ARGV): $_\n";
    }
  }
}

sub check_new_revision {

  if( /^[cC!]\s+(\d{2}\.\d{2}\.\d{4})\s+(\S+)\s+/ ) {
    my $dev = $2;
    my $date = $1;
    my $year = $date;
    $year =~ s/^.*\.//;
    my @devs = split(/\&/,$dev);
    foreach my $d (@devs) {
      #print STDERR "$d   $date  $::file\n" unless $::copy;
      #$::devyear{$d} .= "$year,";
      #$::devcount{$d}++;
    }
    return 1;
  } elsif( /^[cC!]\s+\.\.\.\s+/ ) {	#continuation line
    return 2;
  } else {
    return 0;
  }
}

sub check_old_revision {

  if( /^[cC!].*(\d{2}\.\d{2}\.\d{2})\s+/ ) {
    my $orig = $_;
    if( /^([cC!])\s+(.+)\s+(\d{2}\.\d{2}\.\d{2})\s+(.+)/ or
	/^([cC!])(\s+)(\d{2}\.\d{2}\.\d{2})\s+(.+)/ ) {
      my $comment = $1;
      my $verb = $2;
      my $date = $3;
      my $descrp = $4;
      $verb =~ s/\s+on\s*$//;
      $date = adjust_year($date);
     
      if( $verb eq "written" or $verb eq "revised" or $verb eq "changed" ) {
        ;
      } elsif( $verb =~ /\s*/ ) {	#simply old date
        ;
      } else {
	print STDERR "*** cannot parse (1) ($ARGV): $_\n";
	return 0;
      }
      $descrp =~ s/\s*by ggu\s*//;
      $descrp =~ s/\s*ggu\s*//;
      $_ = "$comment $date\tggu\t$descrp";
      print STDERR "orig: $orig\n";
      print STDERR "new:  $_\n";
      return 1;
    } else {
      print STDERR "*** cannot parse (2) ($ARGV): $_\n";
      return 0;
    }
  } else {
    return 0;
  }
}

sub check_obsolete_revision {

  if( /^[cC!].*(\d{2}\s+\d{2}\s+\d{2,4})/ ) {
    print STDERR "*** obsolete date ($ARGV): $_\n";
    return 1;
  } else {
    return 0;
  }
}

#------------------------------------------------------------------

sub make_date {

  my $date = shift;

  my $year = int($date/10000);
  $date = $date - $year*10000;

  my $month = int($date/100);
  my $day = $date - $month*100;

  $month = "0$month" if $month < 10;
  $day = "0$day" if $day < 10;

  return "$day.$month.$year";
}

sub make_idate {

  my $date = shift;

  if( $date =~ /^(\d\d)\.(\d\d)\.(\d\d\d\d)$/ ) {
    my $idate = $1 + 100*$2 + 10000*$3;
    return $idate
  } else {
    die "cannot parse date ($ARGV): $date\n";
  }
}

sub adjust_year {

  my $date = shift;

  my @f = split(/\./,$date);
  if( $f[2] > 70 ) {
    $f[2] += 1900;
  } else {
    $f[2] += 2000;
  }

  $date = join(".",@f);
  return $date;
}

#------------------------------------------------------------------

sub init_revision {

  return if $::has_been_initialized;

  $::has_copyright = 0;
  $::has_shyfem = 0;
  $::is_manual = 0;
  $::has_been_initialized = 1;
  $::has_revision_log = 0;

  $::comment_char = "";
  $::final_revlog = "";

  %::names = ();
  make_names();
}

#------------------------------------------------------------------

sub make_names {

  %::names = (
         'ggu' => 'Georg Umgiesser'
        ,'aac' => 'Andrea Cucco'
        ,'aar' => 'Aaron Roland'
        ,'ccf' => 'Christian Ferrarin'
        ,'cpb' => ''
        ,'dbf' => 'Debora Bellafiore'
        ,'dmk' => 'Donata Melaku Canu'
        ,'erp' => 'Erik Pascolo'
        ,'fdp' => 'Francesca De Pascalis'
        ,'isa' => 'Isabella Scroccaro'
        ,'ivn' => 'Ivan Federico'
        ,'laa' => 'Leslie Aveytua'
        ,'lcz' => 'Lucia Zampato'
        ,'mbj' => 'Marco Bajo'
        ,'mcg' => 'Michol Ghezzo'
        ,'pzy' => 'Petras Zemlys'
        ,'wmk' => 'William McKiver'
        ,'git' => 'Git Versioning System'
    );
}

#-------------------------------
1;
#-------------------------------

