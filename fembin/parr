#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
# this file has been re-arranged from a version by jgreely@cis.ohio-state.edu
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#parr:
# rearrange conforming PS code to print the pages in an arbitrary
# order.  The -[sS] options (for signature order) assume two-up, left
# to right.  The -o option takes a list of ranges, like this:
#	1-5    1-10,11-20    11-,1-10
# usage: parr [-o list] [-s] [-S n] [file]
#
# jgreely@cis.ohio-state.edu, 89/10/23
# changes by ggu, 96/01/18
# 27.03.2012	ggu	do not split in EPSF section (in_epsf)
# 02.05.2013	ggu	allow for single pages to be extracted
#
# ggu -A to separate pages...

###########################################################################
# set up variables
###########################################################################

$order='';
$signFlag='';
$signCount=0;
$DEBUG=0;
$rangePrint=0;
$TMPDIR='/tmp';

###########################################################################
# determine command line flags
###########################################################################

while ($_ = $ARGV[0],/^-/) {
	shift;
	last if /^-\-$/;
	/^-o/ && ($order = shift,next);
	/^-S/ && ($signCount = shift,$signFlag++,next);
	/^-s/ && ($signFlag++,next);
	/^-d/ && ($DEBUG++,next);
	/^-A/ && ($all=1,next);
	/^-r/ && ($rangePrint++,next);
	die "usage: parr [-A] [-d] [-r] [-o list] [-s] [-S n] [file]\n";
}
if ($signFlag && $order) {
	die "parr: -s and -o cannot be used together\n";
}

###########################################################################
# splitt PS file in "page files" : one for each page + header + trailer
###########################################################################

######## get output file names

$noutfile = @ARGV;
$outfile = $ARGV[0];

#print STDERR "$noutfile     $outfile\n";

if( $noutfile == 1 ) { 		#exactly one file name
  $outfile =~ s/\.ps$//;
  $outfile = "plot" unless $outfile;
} else {			#more file names or stdin
  $outfile = "plot"
}

#print STDERR "$noutfile     $outfile\n";

##############################

$file = "$TMPDIR/p$$.header";
@files = ($file);
$sheet=0;
$in_epsf = 0;
open(file,">$file") ||
  die "$file: $!, stopped";
while (<>) {
	$in_epsf = 1 if /^BeginEPSF/;
	$in_epsf = 0 if /^EndEPSF/;
	#
	# hack to use NeXT Preview: strip old '%%Pages:' lines
	#
	next if /^%%Pages:/;
	if (/^%%Page:/ and not $in_epsf) {
		close(file);
		$sheet++;
		$file = "$TMPDIR/p$$.$sheet";
		push(@files,$file);
		open(file,">$file") ||
		  die "$file: $!, stopped";
	}
	if (/^%%Trailer/ and not $in_epsf) {
		close(file);
		$file = "$TMPDIR/p$$.trailer";
		push(@files,$file);
		open(file,">$file") ||
		  die "$file: $!, stopped";
	}
	print file $_;
}
close(file);

###########################################################################
# determine order in which to print pages
###########################################################################

@order = ();
if ($order) {
	foreach $range (split(/,/,$order)) {
		($start,$sep,$end) = split(/(-)/,$range);
		$start = 1 unless $start;
		$end = $sheet unless $end;
		if ($sep) {
			push(@order,$start..$end);
		}else{
			push(@order,$start);
		}
	}
}elsif ($signFlag) {
	if (! $signCount) {
		$signCount = $sheet;
		$signCount += (4 - $sheet % 4) if ($sheet % 4);
	}else{
		$signCount *=4;
	}
	for($base=0;$base<$sheet;$base+=$signCount) {
		@tmp = ($signCount/2+$base);
		push(@tmp,$tmp[0]+1,$tmp[0]+2,$tmp[0]-1);
		while ($tmp[3] > $base) {
			push(@order,@tmp);
			@tmp = ($tmp[0]-2,$tmp[1]+2,$tmp[2]+2,$tmp[3]-2);
		}
	}
}else{
	@order = (1..$sheet);
}

###########################################################################
# insert blank page if order number is higher than available pages
###########################################################################

@tmp=@order;
@order=();
foreach $page (@tmp) {
	push(@order,$page > $sheet ? "B" : $page);
}

###########################################################################
# do only a range print and exit
###########################################################################

if ($rangePrint) {
	print join(',',@order),"\n";
	unlink @files unless $DEBUG;
	exit(0);
}

unless( $all ) {

###########################################################################
# print header
###########################################################################

open(file,"$TMPDIR/p$$.header");
$_ = <file>;
print $_,"%%Pages: (atend)\n";
print while <file>;
close(file);

###########################################################################
# print pages
###########################################################################

foreach $page (@order) {
	$count++;
	print "%%Page: $count $count\n%%OldPage: $page\n";
	if ($page eq "B") {
		print "showpage\n";
	}else{
		open(file,"$TMPDIR/p$$.$page");
		while (<file>) {
			if( /^%%Page:/ ) {
				s/%%Page:/%%OldPageComment:/;
			}
			print;
		}
		close(file);
	}
}

###########################################################################
# print trailer
###########################################################################

open(file,"$TMPDIR/p$$.trailer");
while (<file>) {
	print unless( /^%%EOF/ );
}
close(file);
print "%%Pages: $count\n%%EOF\n";

} # unless( $all )

###########################################################################
# print single pages
###########################################################################

if( $all ) {
  foreach $page (@order) {
	#print STDERR "$page\n";
	print STDERR " $page";
	open(OUT,">$outfile.$page.ps");

		open(file,"$TMPDIR/p$$.header");
#------------------------------ ggu -> first bounding box
		$_ = <file>;
		print OUT $_;
#------------------------------
		$_ = <file>;
		print OUT $_,"%%Pages: (atend)\n";
		print OUT while <file>;
		close(file);

		print OUT "%%Page: 1 1\n%%OldPage: $page\n";
                open(file,"$TMPDIR/p$$.$page");
                while (<file>) {
                        if( /^%%Page:/ ) {
                                s/%%Page:/%%OldPageComment:/;
                        }
                        print OUT;
                }
                close(file);

		open(file,"$TMPDIR/p$$.trailer");
		while (<file>) {
        		print OUT unless( /^%%EOF/ );
		}
		print OUT "%%Pages: 1\n%%EOF\n";
		close(file);

	close(OUT);
  }
  print STDERR "\n";
}

###########################################################################
# delete tmp files
###########################################################################

unlink @files unless $DEBUG;
exit(0);
