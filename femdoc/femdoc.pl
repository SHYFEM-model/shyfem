#!/usr/bin/perl
#
# extracts documentation from subsys
#
#----------------------------------------------

$indocs = "";

while(<>) {

  if( /^[cC!] DOCS\s+END/ ) { 
	docs_end(); 
  }

  if( $indocs ) {
	docs_process();
  }

  if( /^[cC!] DOCS\s+START\s*(\S*)/ ) { 
	docs_start($1);
  }

}

#----------------------------------------------

sub docs_end {

  if( not $indocs ) {
	print STDERR "No DOCS file opened. Cannot close.\n";
	die "*** Error in structure of documentation.\n";
  }

  $indocs = "";
  close(DOCS);

}

sub docs_start {

  my $file = shift;

  if( $indocs ) {
	print STDERR "Old file $indocs must be closed befor\n";
	print STDERR "new file $file can be opened.\n";
	die "*** Error in structure of documentation.\n";
  }

  $indocs = $file;
  print STDERR "writing file $file.tex\n";
  open(DOCS,">$file.tex");
}

sub docs_process {

  if( /^[cC!] DOCS\s+(\S+)\s+(.+)\n/ ) {	#docs formatting hint
	docs_close_description();
	docs_new_paragraph($2);
  } elsif ( /^cc/ ) {				#normal comment
	docs_close_description();
  } elsif ( /^!!/ ) {				#normal comment
	docs_close_description();
  } elsif ( /^[cC!]\w+#/ ) {			#comment for environments
	docs_close_description();
  } elsif ( /^[cC!](\s*)(.*)\n$/ ) {		#docs comment
	$code = 0;
	docs_format($1,$2);
  } else {					#code or else
	docs_close_description();
	if( $code == 0 ) {
	  print DOCS "\n";
	}
	$code++;
  }
  
}

sub docs_new_paragraph {

  my $paragraph = shift;

  print DOCS "\n\\paragraph\{$paragraph\}\n";
}

sub docs_close_description {

  if( $descrp ) {
	$descrp = "";
	print DOCS "}\n\\par\n";
  }
}

sub docs_open_description {

  my ($descrp,$line) = @_;

  $descrp =~ s/\s+$//;
  $descrp =~ s/ //g;		#remove blanks
  $descrp =~ s/,/\\\\/g;	#substitute comma with new line

  $line =~ s/^\s+//;

  print DOCS "\\descrpitem{$descrp}\n";
  print DOCS "\\descrptext{%\n";
  print DOCS "$line\n";
}

sub docs_continue_description {

  my $line = shift;

  unless( $descrp ) {
	print STDERR "*** Error in line:\n";
	print STDERR "\n";
	print STDERR "$line\n";
	print STDERR "\n";
	print STDERR "This should be a continuation of a description.\n";
	print STDERR "But no description is open.\n";
	print STDERR "Be sure to use tabs for description lines.\n";
	die "*** Error in structure of documentation.\n";
  }
  print DOCS "$line\n";

}

sub docs_format {

  my ($space,$line) = @_;

  if( $space =~ /\t/ ) {		#continue description
	docs_continue_description($line);
  } elsif( $line =~ /(.+)\t(.*)/ ) {	#new description
	docs_close_description();
	$descrp = $1;
	docs_open_description($descrp,$2);
  } else {				#end of description
	docs_close_description();
	print DOCS "$line\n";
  }
}





