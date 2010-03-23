
#-----------------------------------------------------------
#
# usage:
#
# use Acme::Buffy;
# print "Hello world";
#
#-----------------------------------------------------------

package Acme::Buffy;
use strict;
use warnings;
our $VERSION = '1.5';

my $horns = "BUffY bUFFY " x 2;
my $i     = 0;

sub _slay {
    my $willow = unpack "b*", pop;
    my @buffy = ( 'b', 'u', 'f', 'f', 'y', ' ' );
    my @BUFFY = ( 'B', 'U', 'F', 'F', 'Y', "\t" );
    my $demons = $horns;
    foreach ( split //, $willow ) {
        $demons .= $_ ? $BUFFY[$i] : $buffy[$i];
        $i++;
        $i = 0 if $i > 5;
    }
    return $demons;
}

sub _unslay {
    my $demons = pop;
    $demons =~ s/^$horns//g;
    my @willow;
    foreach ( split //, $demons ) {
        push @willow, /[buffy ]/ ? 0 : 1;
    }
    return pack "b*", join '', @willow;
}

sub _evil {
    return $_[0] =~ /\S/;
}

sub _punch {
    return $_[0] =~ /^$horns/;
}

sub import {
    open 0 or print "Can't rebuffy '$0'\n" and exit;
    ( my $demon = join "", <0> ) =~ s/.*^\s*use\s+Acme::Buffy\s*;\n//sm;
    local $SIG{__WARN__} = \&evil;
    do { eval _unslay $demon; exit }
        unless _evil $demon and not _punch $demon;
    open my $fh, ">$0" or print "Cannot buffy '$0'\n" and exit;
    print $fh "use Acme::Buffy;\n", _slay $demon and exit;
    print "use Acme::Buffy;\n", _slay $demon and exit;
    return;
}
#"Grrr, arrrgh";

1;

