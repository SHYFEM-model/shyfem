#!/usr/bin/perl
#
# reformats log of git files
#
#--------------------------------------------------------

make_dev_names();

while(<>) {

  chomp;

  if( /^commit (\w+)/ ) {
    print_commit();
    $commit = $1;
  } elsif( /^Author: (.*)/ ) {
    $author = $1;
    $author =~ s/\s*<.*//;
  } elsif( /^Date:\s+\'(.*)\'/ ) {
    $date = $1;
  } elsif( /^=+ (\S+)/ ) {
    $file = $1;
  } else {
    s/^\s+//; s/\s+$//;
    $comment = $_ if $_;
  }

}

print_commit();		# last commit

print_authors() unless $file;

#--------------------------------------------------------

sub print_commit
{
  return unless $commit;

  $commit = substr($commit,0,10);
  my $acron = make_author_acronym($author);

  if( $file ) {
    my $l = length($file);
    $comment = substr($comment,0,42-1-$l);
    print "$commit $date $acron  $file $comment\n";
  } else {
    $comment = substr($comment,0,42);
    print "$commit $date $acron  $comment\n";
  }
}

sub print_authors
{
  print "Authors of commits:\n";
  for my $author (sort keys %::authors) {
    my $count = $::authors{$author};
    my $fullname = $::dev_names{$author};
    print "  $author \t $count \t $fullname\n";
  }

  return unless scalar %::unknown_authors;

  print "Unknown authors of commits:\n";
  for my $author (sort keys %::unknown_authors) {
    my $count = $::unknown_authors{$author};
    print "  $author  $count\n";
  }
}

#--------------------------------------------------------

sub make_author_acronym
{
  my $author = shift;

  my $acron = $::acronyms{$author};
  unless( $acron ) {
    print STDERR "unknown author: $author\n"; 
    $acron = "nnn";
    $::unknown_authors{$author}++
  }

  $::authors{$acron}++;

  return $acron;
}

#--------------------------------------------------------

sub make_dev_names {

  return if $::dev_initialized;
  $::dev_initialized = 1;

  %::dev_names = (
         'ggu' => 'Georg Umgiesser'
        ,'aac' => 'Andrea Cucco'
        ,'aar' => 'Aaron Roland'
        ,'ccf' => 'Christian Ferrarin'
        ,'cpb' => 'unknown'
        ,'dbf' => 'Debora Bellafiore'
        ,'dmc' => 'Donata Melaku Canu'
        ,'erp' => 'Erik Pascolo'
        ,'fdp' => 'Francesca De Pascalis'
        ,'gir' => 'Ginevra Rosati'
        ,'gml' => 'Giorgio Micaletto'
        ,'clc' => 'Celia Laurent'
        ,'clr' => 'Celia Laurent'
        ,'cla' => 'Carl Amos'
        ,'isa' => 'Isabella Scroccaro'
        ,'ivn' => 'Ivan Federico'
        ,'laa' => 'Leslie Aveytua'
        ,'lcz' => 'Lucia Zampato'
        ,'lrp' => 'Luca Arpaia'
        ,'mbj' => 'Marco Bajo'
        ,'mcg' => 'Michol Ghezzo'
        ,'pzy' => 'Petras Zemlys'
        ,'unm' => 'Urs Neumeier'
        ,'wmk' => 'William McKiver'
        ,'riz' => 'Rasa Idzelyte'
        ,'jca' => 'Jacopo Alessandri'
        ,'git' => 'Git Versioning System'
        ,'shy' => 'The SHYFEM Team'
    );

  %::git_names = (
         'g u' => 'ggu'
        ,'georgu' => 'ggu'
        ,'gmicaletto' => 'gml'
        ,'CeliaLaurent' => 'clr'
        ,'marcobj' => 'mbj'
        ,'christian' => 'ccf'
        ,'dbellafiore' => 'dbf'
        ,'fivan' => 'ivn'
        ,'jalessandri' => 'jca'
        ,'model' => 'shy'
    );

  %::subst_dev_names = (
         'georg' => 'ggu'
        ,'dmk' => 'dmc'
        ,'cl' => 'clc'
        ,'gr' => 'gir'
    );

  foreach my $acron ( keys %::dev_names ) {
    my $name = $::dev_names{$acron};
    $::acronyms{$name} = $acron;
  }
  foreach my $name ( keys %::git_names ) {
    my $acron = $::git_names{$name};
    $::acronyms{$name} = $acron;
  }

}

