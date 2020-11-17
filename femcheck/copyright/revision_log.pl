#!/usr/bin/perl -ws

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

# header
#   header0		empty if no copy, everything in headermain
#   copy		empty lines before and after not counted
#   headermain		everything else
# headermain
#   header1		everything if no revision log
#   rev			start with "revision log :", not excluding this line
#   header2		start with no date line after rev, emtpy if no rev
# body

use lib ("$ENV{SHYFEMDIR}/femcheck/copyright"
		,"$ENV{HOME}/shyfem/femcheck/copyright");

use strict;
use revision_log;
use utils_develop;
use utils_crewrite;

#--------------------------------------------------------------

$::file = $ARGV[0];

my $shydir = "$ENV{SHYFEMDIR}";
my $home = "$ENV{HOME}";
$shydir = "$home/shyfem" unless $shydir;
my $shycopy = "$shydir/femcheck/copyright";

$::copyright_standard = "Copyright (C) 1985-2020  Georg Umgiesser";
$::copyright_full = "$shycopy/copyright_notice.txt";
$::copyright_short = "$shycopy/copyright_short.txt";

$::check = 0 unless $::check;		#checks files
$::gitrev = 0 unless $::gitrev;		#uses gitrev for revision log
$::gitmerge = 0 unless $::gitmerge;	#merges gitrev for revision log
$::stats = 0 unless $::stats;		#uses gitrev for revision log
$::crewrite = 0 unless $::crewrite;	#re-writes c revision log
$::substdev = 0 unless $::substdev;	#substitute developer names
$::onlycopy = 0 unless $::onlycopy;     #only check copyright
$::updatecopy = 0 unless $::updatecopy;	#re-writes copyright section
$::newcopy = 0 unless $::newcopy;       #substitute new copyright

$::copyright = 0;
$::shyfem = 0;
$::manual = 0;
$::copyatend = 0;
$::cstyle_revlog = 0;

$::debug = 0;

#--------------------------------------------------------------

$::type = find_file_type($::file);
$::comchar = get_comment_char($::type);
$::need_revlog = is_code_type($::type);
$::need_revlog = 0 if $::onlycopy;
$::need_copy = is_code_type($::type);
$::need_copy += $::type eq "script";
$::need_copy += $::type eq "text";
$::need_copy += $::type eq "tex";
$::write_file = 1;
$::write_file = 0 if $::onlycopy;

exit 0 if $::type eq "binary";
exit 0 if $::type eq "image";

#print STDERR "reading file $::file with type $::type\n";

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

my ($header,$body) = extract_header();
my ($header0,$copy,$headermain) = extract_copy($header);
my ($header1,$rev,$header2) = extract_revlog($headermain);

if( $::debug ) {
  print_file($header0,$copy,$header1,$rev,$header2);	#use this for debug
}

if( $::crewrite ) {
  ($rev,$header0,$header1,$header2) = 
	handle_c_revlog($rev,$header0,$header1,$header2);
}

check_revlog_dates($rev);
check_copyright($copy);

my $devs = count_developers($rev);
if( ( $::updatecopy or $::newcopy ) and not $::manual ) {
  $copy = handle_copyright($rev,$copy);
}

if( $::substdev and not $::manual ) {
  $rev = subst_developers($rev);
}

if( ( $::gitrev or $::gitmerge ) and ( not $::manual ) ) {
  $header2 = copy_divisor($header1,$copy) unless @$rev;
  if( ( $::gitrev and not @$rev ) or ( $::gitmerge ) ) {
    $rev = integrate_revlog($rev);
  }
}

if( $::stats ) {
  stats_file($::file,$rev,$devs);
} elsif( $::write_file ) {
  write_file("$::file.revnew",$header0,$copy,$header1,$rev,$header2,$body);
}

my $error = $::copyerror;
$error++ if $::need_revlog and not @$rev and not $::manual;
exit $error;

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

sub extract_header
{
  my @header = ();
  my @body = ();

  my $cc = quotemeta($::comchar);
  $cc = "" if $::type eq "special";

  my $in_header = 1;

  while(<>) {

    chomp;

    if( /^\s*$/ ) {			# empty line
    } elsif( $::type eq "fortran" ) {
      if( /^[!cC\*]/ ) {		# comment in first col
      } elsif( /^\s*!/ ) {		# comment with leading blanks
      } else {				# first line of code
        $in_header = 0;
      }
    } elsif( $::type eq "c" ) {
      if( /^\/\*\*/ ) {			# c comment start
      } elsif( /^\\\*\*/ ) {		# c comment end
      } elsif( /^\s*\*/ ) {		# c comment
      } else {				# first line of code
        $in_header = 0;
      }
    } else {
      if( /^\s*$cc/ ) {			# comment
      } else {				# first line not comment
        $in_header = 0;
      }
    }

    if( $in_header ) {
      push(@header,$_);
    } else {
      push(@body,$_);
    }
  }

  return (\@header,\@body);
}

sub extract_copy
{
  my $header = shift;

  my $in_copy = 0;
  my @copy = ();
  my @header0 = ();
  my @headermain = ();

  foreach (@$header) {
    #print "$in_copy: $_\n";
    if( /^..------------------------------/
	or /^\s*.\*\*\*\*\*/
	#or /^\% \*\*\*\*\*/
	) {
      $in_copy++;
      if( $in_copy < 3 ) {
        push(@copy,$_);
        next;
      }
    }
    if( $in_copy == 0 ) {
      push(@header0,$_);
    } elsif( $in_copy == 1 ) {
      $::copyright++ if( /^..\s*Copyright/ );
      $::copyright++ if( /^\% \* Copyright/ );		#for ps files
      $::shyfem++ if( /^..\s*This file is part of SHYFEM./ );
      $::manual++ if( /^..\s*This file is part of SHYFEM.\s+\(m\)/ );
      push(@copy,$_);
    } else {
      push(@headermain,$_);
    }
  }

  if( $::copyright ) {
    if( @header0 > 10 ) {
      if( $::type ne "special" ) {
        print STDERR "    $::file has copyright at end of file\n";
      }
      $::copyatend = 1;
    }
  }

  $::copyerror = 0;

  if( $::copyright ) {
    #if( $::copyright > 1 ) {
    #  print STDERR "*** more than one copyright in file $::file\n";
    #}
    if( not $::stats ) {
      if( $in_copy < 2 ) {
        print STDERR "*** copyright notice not finished in file $::file\n";
        $::copyerror++;
      }
      if( $::shyfem == 0 ) {
        print STDERR "*** missing shyfem line in file $::file\n";
        $::copyerror++;
      }
    }
  } else {
    if( not $::stats and $::need_copy ) {
      print STDERR "*** no copyright in file $::file\n";
      $::copyerror++;
    }
    @headermain = @$header;
    @copy = ();
    @header0 = ();
    if( $::type eq "script" ) {
      push(@header0,shift(@headermain));
    }
  }

  return (\@header0,\@copy,\@headermain);
}

sub extract_revlog
{
  my $headermain = shift;

  my $revs = 0;
  my $in_rev = 0;
  my @rev = ();
  my @header1 = ();
  my @header2 = ();

  foreach (@$headermain) {
    if( /^..\s*revision log :/ ) {
      $in_rev++;
    } elsif( /^..\s*Revision History:/ ) {
      $::cstyle_revlog = 1;
      $in_rev++;
    } elsif( /^[cC!\*]\s*$/ or /^\s*$/ 
		or /^ \*\s*$/ or /^ \*\s+\*\s*$/
		or /^\#\s*$/ ) {
      if( $in_rev == 1 and $revs > 1 ) {
        $in_rev++;				#account for empty line
      }
    }
    if( $in_rev == 0 ) {
      push(@header1,$_);
    } elsif( $in_rev == 1 ) {
      $revs++;
      push(@rev,$_);
    } else {
      push(@header2,$_);
    }
  }

  if( $in_rev > 2 ) {
    print STDERR "*** more than one revision log in file $::file\n";
  }
  if( $::cstyle_revlog ) {
    print STDERR "*** old c style revision log in file $::file\n";
  }
  if( $::need_revlog and not $revs and not $::manual ) {
    print STDERR "*** no revision log in file $::file\n";
  }

  return (\@header1,\@rev,\@header2);
}

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

sub read_file
{
  my $file = shift;

  open(FILE,"<$file") || die "*** cannot open file $file\n";;
  my @text = <FILE>;
  close(FILE);

  return \@text;
}

sub write_file
{
  my $file = shift;

  open(FILE,">$file");

  foreach my $block (@_) {
    foreach my $line (@$block) {
      print FILE "$line\n";
    }
  }

  close(FILE);
}

sub print_file
{
  my $n = 0;

  my @title = ("header0","copy","header1","rev","header2");

  foreach my $block (@_) {
    $n++;
    my $title = shift(@title);
    print "======================================================\n";
    print "$n - $title\n";
    print "======================================================\n";
    foreach my $line (@$block) {
      print "$line\n";
    }
  }
}

sub stats_file
{
  my ($file,$rev,$dev) =@_;

  my $revs = @$rev;
  $revs -= 2 if $revs > 0;

  my $error = 0;
  $error++ if $::copyright != 1;
  $error++ if $::shyfem != 1;
  $error++ if $::manual != 0;
  $error++ if $revs == 0;

  return unless $error;

  print" $::copyright $::shyfem $::manual $revs $devs   $file\n";
}

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

sub handle_copyright
{
  my ($rev,$copy) = @_;

  if( check_copyright($copy) ) {	# this is true if copy and unknown type
    print "    not changing copyright...\n";
    return $copy;
  }

  if( $::newcopy and $::copyatend ) {
    print "*** $::file has copyright at end of file... use --updatecopy\n";
    return $copy;
  }

  if( $::updatecopy ) {
    return $copy unless @$copy;			#do not update if no copyright
  } elsif( $::newcopy ) {
    $copy = [];					#make a new copyright
  }

  if( not @$copy ) {				# no copyright - integrate
    $copy = integrate_copyright();
  }

  if( @$rev ) {					# has revision log
    my $newcopy = handle_developers($rev,$copy);
    $copy = substitute_copyright($copy,$newcopy);
  } else {					# no revision log - use standard
    my @newcopy = ();
    my $line = $::comchar . "    " . $::copyright_standard;
    push(@newcopy,$line);
    $copy = substitute_copyright($copy,\@newcopy);
  }

  $::copyerror = 0;				# when we get here no error
  return $copy;
}

sub integrate_copyright
{
  my ($copy,$divisor_start,$divisor_end);

  $divisor_start = $::comchar . "-" x 74;
  $divisor_end   = $::comchar . "-" x 74;
  if( $::type eq "c" ) {
    $divisor_start = "/"  . "*" x 72 . "\\";
    $divisor_end   = "\\" . "*" x 72 . "/";
  }

  if( $::type eq "fortran" or $::type eq "c" ) {
    $copy = read_file($::copyright_full);
  } else {
    $copy = read_file($::copyright_short);
  }

  foreach (@$copy) {
    chomp;
    $_ = $::comchar . $_;
  }
  unshift(@$copy,$divisor_start);
  unshift(@$copy,"");
  push(@$copy,$divisor_end);

  return $copy;
}

sub substitute_copyright
{
  my ($cold,$csubst) = @_;

  my @cnew = ();
  my $in_old = 0;

  foreach (@$cold) {
    if( /..\s*Copyright/ ) {
      push(@cnew,@$csubst) unless $in_old;
      $in_old++;
      next;
    }
    push(@cnew,$_);
  }

  return \@cnew;
}

sub check_copyright
{
  my ($copy) = @_;

  return 0 if $::manual;

  if( @$copy and $::type eq "unknown" ) {
    print "*** $::file has copyright but has unknown type\n";
    return 1;
  }
  return 0;
}

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

sub integrate_revlog
{
  my $rev = shift;

  my $nrev = @$rev;
 
  if( $nrev ) {
    print STDERR "*** not ready for merging, can only integrate in $::file\n";
    return $rev;
  }

  my @gitrev=`git-file -revlog $::file`;

  my @new = ();
  my $in_revision_log = 0;
  foreach (@gitrev) {
    $in_revision_log = 1 if /revision log/;
    next unless $in_revision_log;
    chomp;
    push(@new,$_);
  }

  my $nnew = @new;
  if( $nnew <= 3 ) {
    print STDERR "*** no git revlog found in $::file ($nnew)\n";
    return $rev;
  } else {
    print STDERR "    from git $nnew lines read for file $::file\n";
    if( $::type eq "c" ) {
      @new = revadjust_for_c(@new);
    }
  }

  return \@new;
}

sub check_revlog_dates
{
  my $rev = shift;

  if( $::cstyle_revlog ) {
    #print "   cannot check dates in old c style revlog\n";
    return
  }

  my $old_idate = 0;
  foreach (@$rev) {
    my ($date,$dev,$text) = parse_revision_line($_);
    if( $date ) {
      my $idate = make_idate($date);
      if( $old_idate > $idate ) {
        print STDERR "   *** error in dates of revision log: "
				. "$old_idate - $idate ($::file)\n";
      }
      $old_idate = $idate;
    } else {
      if( $dev eq "error" ) {
        print STDERR "   *** error parsing revision log: $_ ($::file)\n";
      }
    }
  }
}

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

sub copy_divisor
{
  my ($header,$copy) = @_;

  my @aux = reverse(@$header);
  my @new = ();
  my $div;

  foreach (@aux) {
    $div = $_;
    last if /^[cC!]\s*\*\*\*/;
    last if /^[cC!]\s*---/;
    last if /^[cC!]\s*===/;
    last if /^\s*\*\*\*/;
    $div = "";
  }

  if( $div ) {
    push(@new,$div);
    push(@new,"\n");
  } else {
    $::divisor_start = $copy->[0];
    $::divisor_end = $copy->[-1];
  }

  return \@new;
}

sub write_array
{
  my $text = shift;

  print "$text\n";
  foreach (@_) {
    print "$_\n";
  }
}

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

