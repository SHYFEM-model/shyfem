#!/usr/bin/perl
#
# extracts documentation from subsys

$indocs = "";

while(<>) {

  if( /^c DOCS\s+END/ ) { &docsend; }

  if( $indocs ) {
	&docsprocess;
  }

  if( /^c DOCS\s+START\s*(\S*)/ ) { 
	$file = $1;
	&docsstart;
  }

}



#########################################

sub docsend {

  $indocs = "";
  close(DOCS);

}

sub docsstart {

  if( $indocs ) {
	print "Old file $indocs must be closed befor\n";
	print "new file $file can be opened.\n";
	die "Error in structure of documentation.\n";
  }
  $indocs = $file;
  $section = "";
  print STDERR "writing file $indocs.tex\n";
  open(DOCS,">$indocs.tex");
}

sub docsprocess {

  if( /^c DOCS\s+(\S+)\s+(.+)\n/ ) {	#docs formatting hint
	$section = $1;
	$paragraph = $2;
	&docsclose;
	&docspara;
  } elsif ( /^cc/ ) {			#normal comment
	&docsclose;
  } elsif ( /^c\w+#/ ) {		#comment for environments
	&docsclose;
  } elsif ( /^c(\s*)(.*)\n$/ ) {	#docs comment
	$space = $1;
	$line = $2;
	$code = 0;
	&docsformat;
  } else {				#code or else
	&docsclose;
	if( $code == 0 ) {
	  print DOCS "\n";
	}
	$code++;
  }
  
}

sub docspara {

  print DOCS "\n\\paragraph\{$paragraph\}\n";
}

sub docsclose {

  if( $descrp ) {
	$descrp = "";
	print DOCS "}\n\\par\n";
  }
}

sub docsformat {

  if( $space =~ /\t/ ) {		#tab as first space
	die "internal error (1): $line" unless $descrp;
  	print DOCS "$line\n";
  } elsif( $line =~ /(.+)\t(.*)/ ) {	#new description
	&docsclose;
	$descrp = $1;
	$line = $2;
	$descrp =~ s/\s+$//;
	$descrp =~ s/ //g;		#remove blanks
	$descrp =~ s/,/\\\\/g;		#substitute comma with new line
	$line =~ s/^\s+//;
	print DOCS "\\descrpitem{$descrp}\n";
	print DOCS "\\descrptext{%\n";
	print DOCS "$line\n";
  } else {
	&docsclose;
	print DOCS "$line\n";
  }
}

