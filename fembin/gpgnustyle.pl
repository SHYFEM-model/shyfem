#!/usr/bin/perl
#
# changes style in gnu plot

$type = shift;
@psfile = <>;

#---------------------------------- first pass

foreach (@psfile) {

  if( /^(LT\d)/ ) {
    &close_old;
    &open_new;			#here data is set it is the number of LT
    $special{$data} = 0;
  } elsif( /^grestore/ ) {
    &close_old;
  }

  if( $data ) {
    if( / Box/ ) {		#this is how we identify special plot
      $special{$data} = 1;
    }
  }

}
    
foreach (keys %special) {
  $howmany += $special{$_};
}
if( $howmany == 0 ) {	#set last style encountered
  $special{$lastdata} = 1;
  print STDERR "last data set to color... $lastdata\n";
}

#---------------------------------- second pass

@new = ();
$writing = 1;

foreach (@psfile) {

  if( /^(LT\d)/ ) {
    &close_old;
    &open_new;
  } elsif( /^grestore/ ) {
    &close_old;
  }

  if( $data && $special{$data} ) {
    if( /^$data/ ) {
      $_ = "LT7\n";
    } elsif( /Rshow$/ ) {
      $_ .= "/Color true def LT7\n";
    } elsif( / Box/ ) {
      s/ Box/ BoxF/;
    }
  }

  push(@new,$_);

}
    
print @new;

#---------------------------------- subroutines

sub close_old {

  if( $data ) {
    #print STDERR "closing old data set $data ($special{$data})\n";
    if( $special{$data} && $writing ) {
      push(@new,"/Color false def\n");
    }
    $data = "";
  }
}

sub open_new {

  if( $data ) {
    #print STDERR "old data set not closed: $data\n";
  } else {
    $data = $1;
    $lastdata = $data;
    #print STDERR "opening new data set $data\n";
  }
}

