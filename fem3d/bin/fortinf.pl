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
# analyses fortran 77 files and writes statistics
#
# still to do: 
#
# - handle continuation lines
# - write also dimension of arrays
# - programm
# - * is a comment - ok
# - block data - ok

$ftypes = " logical integer real double doubleprecision complex character ";
$file = "";
$proc = "";
$actproc = "";

while(<>) {

  if ( $file ne $ARGV ) {
    $file = $ARGV;
    print "-------------------------------------\n";
    print "file: $file\n";
    print "-------------------------------------\n";
    print STDERR "$file\n";
  }

  chomp;
  tr/A-Z/a-z/;		#make lowercase
  next if /^c/i;	#traditional comment
  next if /^\*/i;	#traditional comment (unusual)
  if( /^(.*?)!(.*)/ ) {		#new comment till end of line
	$code = $1;
	$rest = $2;
	unless( $rest =~ /\'/ ) {
		$_ = $code . "\n";
	}
  }

  &untab;

  $c = substr($_,5,1);
  if( $c && $c ne " " ) {
    $contilines{$file}++;
    #print STDERR "conti: $_\n";
  }

  $found = "";

  s/double\s+precision/doubleprecision/;
  s/block\s+data/blockdata/;

  if( $actproc ) {
    if( /^\s+end\s*$/ ) {
	$found = "end";
    } else {
	$found = "";
    }
  } else {
    if( /^\s+program\s+(\w*)/ ) {
	$found = "program";
	$actual = $found;
	$proc = $1;
    } elsif( /^\s+subroutine\s+(\w*)/ ) {
	$found = "subroutine";
	$actual = $found;
	$proc = $1;
    } elsif( /^\s+function\s+(\w*)/ ) {
	$found = "function";
	$actual = $found;
	$proc = $1;
    } elsif( /^\s+(\w+)\s+function\s+(\w*)/ ) {
      $type = $1;
      $second = $2;
      if( $ftypes =~ /\s+$type\s+/ ) {
	$found = "function";
	$actual = $found;
	$proc = $second;
        $typedfunction{$file} .= "$proc ";	#function typed on same line
      } else {
	die "bogus function declaration: $_\n";
	$found = "";
      }
    } elsif( /^\s+blockdata\s+(\w*)/ ) {
	$found = "blockdata";
	$actual = $found;
	$proc = $1;
    } elsif( /^\s*$/ ) { 	#empty line -> ok
	$found = "";
    } else {
	print STDERR "not in procedure: $file -> main\n$_\n";
	die "not continuing...\n";
	$found = "program";
	$actual = $found;
	$proc = "main";
	push(@noproc,$file);
    }
  }

  if( $found ) {	#procedure declaration or similar
    if( $found eq "end" ) {
	unless( $actproc ) {
	  #die "not in procedure: $file\n$_";
	  print STDERR "not in procedure: $file -> main\n$_";
	  push(@noproc,$file);
	  $proc = "main";
	  die "this should never happen...\n";
	}
	push(@proc,"$_\n");
	processproc();		#here elaborate procedure
	$actual = "";
	$proc = "";
	$actproc = "";
    } else {		#procedure declaration
	die "already in procedure: $file $proc $actproc\n$_" if $actproc;
	$actproc = $proc;
	@proc = ();
	push(@proc,$_);
    }
  } elsif( $actproc ) {	#we are in a procedure -> ok
	# nothing
  } else {		#we are outside a procedure
    	if( /\S+/ ) {
	    print STDERR "statement outside procedure... \"$_\"\n";
    	}
  }
  if( $actual ) {
	push(@proc,$_);
  }

}

print "\n\n\n";

make_hash_unique( \%commonfile );
make_hash_unique( \%commonproc );
make_hash_unique( \%impliproc );
make_hash_unique( \%typedfunction );

open(COMMON,">COMMONF");
print "=======================================================\n";
print "commonblocks by file:\n";
print "=======================================================\n";
foreach $common (keys %commonfile) {
  $line = $commonfile{$common};
  @common = &get_values($line);
  $n = @common;
  print "/$common/ ($n) $line\n";
  print COMMON "/$common/ ($n) $line\n";
}
close(COMMON);

print "=======================================================\n";
print "commonblocks by procedure:\n";
print "=======================================================\n";
foreach $common (keys %commonproc) {
  $line = $commonproc{$common};
  @common = &get_values($line);
  $n = @common;
  print "/$common/ ($n) $line\n";
}

print "=======================================================\n";
print "implicit declarations:\n";
print "=======================================================\n";
foreach $file (keys %impliproc) {
  $line = $impliproc{$file};
  @key = &get_values($line);
  $n = @key;
  print "$file ($n) $line\n";
}

print "=======================================================\n";
print "typed functions:\n";
print "=======================================================\n";
foreach $file (keys %typedfunction) {
  $line = $typedfunction{$file};
  @key = &get_values($line);
  $n = @key;
  print "$file ($n) $line\n";
}

print "=======================================================\n";
print "continuation lines:\n";
print "=======================================================\n";
foreach $file (keys %contilines) {
  print "$file: $contilines{$file}\n";
}

print "=======================================================\n";
print "no procedure declaration:\n";
print "=======================================================\n";
foreach $file (@noproc) {
  print "$file\n";
}

check_modules( \%commonproc );

#####################################################################

sub make_unique {

  my $line = shift;

  my @f = split(" ",$line);
  my %f = ();

  foreach my $elem (@f) {
    $f{$elem} = 1;
  }

  @f = ();
  foreach my $elem (keys %f) {
    push(@f,$elem);
  }

  return @f;
}

sub get_values {

  my $line = shift;

  my @f = split(" ",$line);

  return @f;
}

sub make_hash_unique {		#makes values in hash table unique

  my $rhash = shift;

  foreach $key (keys %$rhash) {
    my $value = $$rhash{$key};
    my @f = &make_unique($value);
    my $line = join(" ",@f);
    $$rhash{$key} = $line;
  }
}

#####################################################################

sub print_line {

  $first = shift;
  $_ = shift;
  $max = 60;

  @f = split;

  $aux = $first;
  while( $next = shift(@f) ) {
    $new = "$aux $next";
    $n = length($new);
    #print "$n $max: $new\n";
    if( $n > $max ) {
      if( $aux eq $first ) {
        $old = $new;
	$aux = $first;
      } else {
        $old = $aux;
	$aux = "$first $next";
      }
      print "$old\n";
    } else {
      $aux = $new;
    } #if
  } #while

  if( $aux ne $first ) {
      print "$aux\n";
  }
}

#####################################################################

sub processproc {

  print "=======================================================\n";
  print "$actual $proc\n";
  print "=======================================================\n";

  $max = @proc;
  $impli = 0;
  %type = ();
  %bound = ();
  %funct = ();
  %common = ();

# loop on routine

  for($i=0;$i<$max;$i++) {

    $_ = $proc[$i];	#next line
    s/^\s+//;

    # get first word of line

    if( /^(\w+)\s+/ ) {
	$first = $1;
    } elsif( /^(\w+)\*\S+\s+/ ) {
	$first = $1;
    } else {
	next;	#nothing found -> empty line or else
    }

    if( $ftypes =~ /\s+$first\s+/ ) {
	&type_declaration($first,$_);
    } elsif( $first eq "implicit" ) {
	&implicit($first,$_);
    } elsif( $first eq "common" ) {
	&common($first,$_);
    }
  }

# print type declarations of fucntion (if applicable)

  if( $actual eq "function" ) {
    $funct{$proc} = ":function";
    $var = $proc;
    print "$type{$var} $var$bound{$var}  $funct{$var}\n";
    $type{$var} = "";
  }

# print type declarations

  %types = ();
  foreach $var (keys %type) {
    $type = $type{$var};
    next unless $type;
    $types{$type} .= "$var$bound{var} ";
  }

  foreach $type (keys %types) {
    $line = $types{$type};
    #print "$type $line\n";
    &print_line($type,$line);
  }

# print common blocks

  if( %common ) {
    $line = "";
    foreach $var (keys %common) {
      $line .= "/$var/ ";
    }
    #print "common $line\n";
    &print_line("common",$line);
  }

# memorize if no implicit none

  unless( $impli ) {
    $impliproc{$file} .= "$proc ";
  }

}
  
#####################################################################

sub type_declaration {

  $first = shift;
  $_ = shift;

	#print "$first - $_";

	while( /\((\w+,)+\w+\)/ ) {	#set tempor. (,) to (;)
	  $f = $`;
	  $l = $';
	  $par = $&;
	  $par =~ s/,/;/g;
	  $_ = $f . $par . $l;
	}
	s/,/ /g;			#now delete remaining ,
	s/;/,/g;			#and restore , from ;
	s/^\s*//;			#delete leading white space
	@f = split(" ",$_);
	$first = shift(@f);		#shift off first element (type)
	shift(@f) if $first eq "double";#next word is precision
	
	foreach $i (@f) {
	  if( $i =~ /^\w+/ ) {
	    $var = $&;
	    $bound = $';
	  }
	  $type{$var} = $first;
	  $bound{$var} = $bound;
	  $funct{$var} = "";
	}
  }

sub implicit {

  $first = shift;
  $_ = shift;

    if( /^implicit\s+none\s*\n$/ ) {
	print "implicit none\n";
	$impli = 1;
    }
}

sub common {

  $first = shift;
  $_ = shift;

  s/common\s*//;

  #print "common block: $_";

  if( /^[^\/]+/ ) {	#if first char is not / => blank common
    $name = " ";
    $_ = $';
    $common{$name} = 1;
  }

  while( /^\/\s*(\w+)\s*\/[^\/]*/ ) {	# / name / more
    $name = $1;
    $_ = $';
    $common{$name} = 1;
    $commonfile{$name} .= "$file ";
    $commonproc{$name} .= "$proc ";
    #print "...$name\n";
  }

}

###############################################################

sub untab {

  $tablength = 8;

  while( ($i = index($_,"\t")) != -1 ) {
    $j = $i % $tablength;
    $j = $tablength - $j;
    if( $j > 0 ) {
	substr($_,$i,1) = substr('           ',0,$j);
    }
  }
}
  
###############################################################

sub check_modules {

  my $hash = shift;

  my %commons = ();
  my %commons_left = ();

  my %modules =	(
		 "basin"	=>	"nen3v xgv ygv hm3v hev hkv \
						ipv ipev iarv \
						nkonst descrr"
		,"hydro"	=>	"unv vnv uov vov zov znv zeov zenv \
						ulov vlov ulnv vlnv \
						utlov vtlov utlnv vtlnv \
						wlov wlnv \
						uprv vprv wprv \
						xv up0v vp0v"
		,"dry"		=>	"iwegv"
		,"femtime"	=>	"femtim"
		,"simul"	=>	"descrp"
		,"vertical"	=>	"ilhv ilhkv level hlv hldv"
		,"femmath"	=>	"mkonst"
		,"femphys"	=>	"pkonst fcorv"
		,"links"	=>	"lenkv linkv ilinkv"
		,"aux"		=>	"v1v v2v v3v ve1v \
						saux1 saux2 saux3 saux4"
		,"geom"		=>	"ev ieltv kantv dxv dyv rdistv"
		,"baroc"	=>	"saltv tempv rhov bpresv rtv rsv"
		,"bound"	=>	"bnd irv rqv rzv inodv \
						boundn saltn tempn conzn \
						bio2dn \
						ruv rvv crad \
						rhv rlv rrv ierv \
						iopbnd"
		,"diffusion"	=>	"difv visv difhv austv wdifhv"
		,"depth"	=>	"hdenv hdeov hdknv hdkov hlhv areakv"
		,"conz"		=>	"cnv rcv"
		,"flux"		=>	"kfluxc iflux"
		,"volume"	=>	"kvolc ivol"
		,"wind"		=>	"tauxnv tauynv wxnv wynv wxov wyov \
						ppv pov pnv \
						winpar wintim windat"
		,"tide"		=>	"xgeov ygeov zeqv tidcom astro"
		,"wave"		=>	"waveh wavep waved"
		,"friction"	=>	"czv"
		,"param"	=>	"parsco fnmsco secsec nrdcom"
		,"regular"	=>	"ppp20 mercom mcoo"
		,"ousfile"	=>	"ouscom ousvar"
		,"extfile"	=>	"extcom knausc"
		,"nosfile"	=>	"noscom nosvar"
		,"debug"	=>	"comdebug"
		,"implicit"	=>	"semimi semimr istotc"
		,"grid"		=>	"vscale"
		,"date"		=>	"dtsdts"
		,"histo"	=>	"ihisto rhisto"
		,"gotm"		=>	"nuhv numv tken eps rls"
		,"list"		=>	"lstloc"
		,"light"	=>	"luxlen_c"
		,"system"	=>	"amat"
		,"change"	=>	"art"
		);

# dxv,dyv is not needed
# 
#
#		,"lagrange"	=>	"nbdy est xst yst zst ie_body \
#						x_body y_body z_body tin \
#						tdecay fx vel_ie \
#						lt_body dvert fall alt"

  my $imod = 0;
  my $icom = 0;
  foreach my $line (values %modules) {
    $imod++;
    my @f = split(/\s+/,$line);
    foreach my $common (@f) {
      $icom++;
      $commons{$common} = 1;
    }
  }

  foreach $common (keys %$hash) {
    if( $commons{$common} ) {
      $commons{$common}++;
    } else {
      $commons_left{$common} = $hash->{$common};
    }
  }

  my $ictot = 0;
  foreach $key (keys %commons) {
    my $count = $commons{$key} - 1;
    my $line = "";
    if ( $count <= 0 ) {
      $ictot++;
      $line = "******";
    }
    print STDERR "modules... $key    $line\n";
  }
  print STDERR "number of modules: $imod    number of common handled: $icom\n";
  print STDERR "number of common blocks not referenced: $ictot\n";

  @by_count = ();

  my $ntot = 0;
  my $ctot = 0;
  my $nmax = 0;
  open(COMMON,">COMMONL");
  foreach $common (keys %commons_left) {
    @common = &make_unique($commons_left{$common});
    $n = @common;
    $by_count[$n] .= "$common ";
    $ntot += $n;
    $nmax = $n if $n > $nmax;
    $ctot++;
    $line = join(" ",@common);
    print COMMON "/$common/ ($n) $line\n";
  }
  close(COMMON);

  for( $i=0; $i<=$nmax ; $i++ ) {
    my $com = $by_count[$i];
    if( $com ) {
      print STDERR "$i : $com\n";
    }
  }
    
  print STDERR "common blocks left:   $ctot   $ntot\n";
}

