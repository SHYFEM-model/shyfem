## package Text::PDF::API::Color;

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

## package Graphics::ColorObject;

package Color::Object;

$VERSION='0.1_02';

use POSIX;

=pod

=head1 NAME

Color::Object - A OO-Color Module 

=head1 DESCRIPTION

A module for manipulation Colors within RGB, HSV and HSL color-spaces for
usage within PDF-Documents especially with the Text::PDF::API modules.

=head1 SYNOPSIS

	use Color::Object;

	$cl = Color::Object->new;
	$cl = Color::Object->newRGB($r,$g,$b);
	$cl = Color::Object->newHSV($h,$s,$v);
	$cl = Color::Object->newHSL($h,$s,$l);

	$cl->setRGB($r,$g,$b);
	$cl->addBrightness($br);
	($h,$s,$l) = $cl->asHSL;

=head1 METHODS

=cut

sub mMin {
	my $n=10000000000000;
	map { $n=($n>$_) ? $_ : $n } @_;
	return($n);	
}

sub mMax {
	my $n=0;
	map { $n=($n<$_) ? $_ : $n } @_;
	return($n);	
}

sub RGBtoHSV ($$$) {
	my ($r,$g,$b)=@_;
	my ($h,$s,$v,$min,$max,$delta);

	$min= mMin($r,$g,$b);
	$max= mMax($r,$g,$b);

        $v = $max;                              

        $delta = $max - $min;

        if( $delta > 0.000000001 ) {
                $s = $delta / $max;
        } else {
                $s = 0;
                $h = 0;
                return($h,$s,$v);
        }

        if( $r == $max ) {
                $h = ( $g - $b ) / $delta; 
        } elsif( $g == $max ) {
                $h = 2 + ( $b - $r ) / $delta; 
        } else {
                $h = 4 + ( $r - $g ) / $delta;
	}
        $h *= 60;
        if( $h < 0 ) {$h += 360;}
        return($h,$s,$v);
}
sub HSVtoRGB ($$$) {
	my ($h,$s,$v)=@_;
	my ($r,$g,$b,$i,$f,$p,$q,$t);

        if( $s == 0 ) {
                ## achromatic (grey)
                return ($v,$v,$v);
        }

        $h /= 60;                       ## sector 0 to 5
        $i = POSIX::floor( $h );
        $f = $h - $i;                   ## factorial part of h
        $p = $v * ( 1 - $s );
        $q = $v * ( 1 - $s * $f );
        $t = $v * ( 1 - $s * ( 1 - $f ) );

	if($i<1) {
		$r = $v;
                $g = $t;
                $b = $p;
	} elsif($i<2){
		$r = $q;
                $g = $v;
                $b = $p;
	} elsif($i<3){
		$r = $p;
                $g = $v;
                $b = $t;
	} elsif($i<4){
		$r = $p;
                $g = $q;
                $b = $v;
	} elsif($i<5){
		$r = $t;
                $g = $p;
                $b = $v;
	} else {
		$r = $v;
                $g = $p;
                $b = $q;
	}
	return ($r,$g,$b);
}
sub RGBtoHSL ($$$) {
	my ($r,$g,$b)=@_;
	my ($h,$s,$v,$l,$min,$max,$delta);

	$min= mMin($r,$g,$b);
	$max= mMax($r,$g,$b);
	($h,$s,$v)=RGBtoHSV($r,$g,$b);
	$l=($max+$min)/2.0;
        $delta = $max - $min;
	if($delta<0.00000000001){
		return(0,0,$l);
	} else {
		if($l<=0.5){
			$s=$delta/($max+$min);
		} else {
			$s=$delta/(2-$max-$min);
		}
	}
	return($h,$s,$l);
}
sub RGBquant ($$$) {
	my($q1,$q2,$h)=@_;
	while($h<0){$h+=360;}
	$h%=360;
	if ($h<60) {
		return($q1+(($q2-$q1)*$h/60));
	} elsif ($h<180) {
		return($q2);
	} elsif ($h<240) {
		return($q1+(($q2-$q1)*(240-$h)/60));
	} else {
		return($q1);
	}
}
sub HSLtoRGB ($$$) {
	my($h,$l,$s,$r,$g,$b,$p1,$p2)=@_;
	if($l<=0.5){
		$p2=$l*(1+$s);
	} else {
		$p2=$l+$s-($l*$s);
	}
	$p1=2*$l-$p2;
	if($s<0.0000000000001){
		$r=$l; $g=$l; $b=$l;
	} else {
		$r=RGBquant($p1,$p2,$h+120);
		$g=RGBquant($p1,$p2,$h);
		$b=RGBquant($p1,$p2,$h-120);
	}
	return($r,$g,$b);
}

=item Color::Object->new

=cut 

sub new {
	my $class=shift @_;
	my $self={};
	bless($self,$class);
	return($self);
}

=item Color::Object->newRGB $r, $g, $b

=cut 

sub newRGB {
	my $class=shift @_;
	my ($r,$g,$b)=@_;
	my $self=$class->new;
	$self->setRGB($r,$g,$b);
	return $self;	
}

=item Color::Object->newHSV $h, $s, $v

=cut 

sub newHSV {
	my $class=shift @_;
	my ($h,$s,$v)=@_;
	my $self=$class->new;
	$self->setHSV($h,$s,$v);	
	return $self;	
}

=item Color::Object->newHSL $h, $s, $l

=cut 

sub newHSL {
	my $class=shift @_;
	my ($h,$s,$l)=@_;
	my $self=$class->new;
	$self->setHSL($h,$s,$l);	
	return $self;	
}

=item Color::Object->newGrey $grey

=cut 

sub newGrey {
	my $class=shift @_;
	my ($g)=@_;
	my $self=$class->new;
	$self->setGrey($g);	
	return $self;	
}

=item ( $r, $g, $b ) = $cl->asRGB 

Returns $cl's rgb values. Range [0 .. 1].

=cut 

sub asRGB {
	my $self=shift @_;
	return @{$self->{'rgb'}};
}

=item ( $h, $s, $v ) = $cl->asHSV 

Returns $cl's hsv values. Ranges h [0 .. 360], s/v [0 .. 1].

=cut 

sub asHSV {
	my $self=shift @_;
	return @{$self->{'hsv'}};
}

=item ( $h, $s, $l ) = $cl->asHSL 

Returns $cl's hsl values. Ranges h [0 .. 360], s/l [0 .. 1].

=cut 

sub asHSL {
	my $self=shift @_;
	return @{$self->{'hsl'}};
}

=item $grey = $cl->asGrey 

=item $grey = $cl->asGrey2

Returns $cl's grey value. Range [0 .. 1]. Functions 2 returns the geometric mean of the corresponding RGB values.

=cut 

sub asGrey {
	my $self=shift @_;
	return $self->{'grey'};
}

sub asGrey2 {
	my $self=shift @_;
	my ($r,$g,$b)=@{$self->{'rgb'}};
	return((($r**2+$g**2+$b**2)**0.5)/3);
}

=item ( $c, $m, $y )= $cl->asCMY 

Returns $cl's cmy values. Range [0 .. 1]. 

=cut 

sub asCMY {
	my $self=shift @_;
	return(map { 1-$_ } $self->asRGB);
}

=item ( $c, $m, $y, $k )= $cl->asCMYK 

=item ( $c, $m, $y, $k )= $cl->asCMYK2

=item ( $c, $m, $y, $k )= $cl->asCMYK3

Returns $cl's cmyk values. Range [0 .. 1]. 
Function 2 returns a 25% lighter color-equivalent.
Function 3 returns a 25% lighter color-equivalent.

=cut 

sub asCMYK {
	my $self=shift @_;
	my @cmy=(map { 1-$_ } $self->asRGB);
	my $k=mMin(@cmy);
	return((map { $_-$k } @cmy),$k);
}

sub asCMYK2 {
	my $self=shift @_;
	my @cmyk=$self->asCMYK;
	$cmyk[3]*=0.75;
	return(@cmyk);
}

sub asCMYK3 {
	my $self=shift @_;
	my @cmyk=$self->asCMY;
	$cmyk[3]=0;
	return(map { $_*0.75 } @cmyk);
	return(@cmyk);
}

=item $hex = $cl->asHex 

Returns $cl's rgb values as 6 hex-digits.

=cut 

sub asHex {
	my $self=shift @_;
	return sprintf('%02X%02X%02X',map {$_*255} $self->asRGB);
}

=item $cl->setRGB $r, $g, $b 

Sets the $cl's rgb values. Valid range [0 .. 1].

=cut 

sub setRGB {
	my $self=shift @_;
	my ($r,$g,$b)=@_;
	$self->{'rgb'}=[$r,$g,$b];	
	$self->{'hsv'}=[RGBtoHSV($r,$g,$b)];	
	$self->{'grey'}=(0.299*$r)+(0.587*$g)+(0.144*$b);
	$self->{'hsl'}=[RGBtoHSL($r,$g,$b)];	
}

=item $cl->setHSV $h, $s, $v 

Sets the $cl's hsv values. Valid ranges: h [0..360], s/v [0..1]. 

=cut 

sub setHSV {
	my $self=shift @_;
	my ($h,$s,$v)=@_;
	$self->setRGB(HSVtoRGB($h,$s,$v));	
}

=item $cl->setHSL $h, $s, $l 

Sets the $cl's hsl values. Valid ranges: h [0..360], s/l [0..1]. 

=cut 

sub setHSL {
	my $self=shift @_;
	my ($h,$s,$l)=@_;
	$self->setRGB(HSLtoRGB($h,$s,$l));	
}

=item $cl->setGrey $grey 

Sets the $cl's grey value. Valid range [0 .. 1].

=cut 

sub setGrey {
	my $self=shift @_;
	my ($g)=@_;
	$self->setRGB($g,$g,$g);
}

=item $cl->setHex $hex 

Sets the $cl's rgb values using 6 hex-nibbles.

=cut 

sub setHex {
	my $self=shift @_;
	my ($hx)=@_;
	my($r,$g,$b) = map { $_/255 } unpack('H3',$hx);
	$self->setRGB($r,$g,$b);	
}

=item $cl->addSaturation $saturation 

Adds to the $cl's saturation in the HSV model. Valid range [-1 .. 1].

=cut 

sub addSaturation {
	my $this=shift @_;
	my $sat=shift @_;
	my ($h,$s,$v)=$this->asHSV;
	$this->setHSV($h,$s+$sat,$v);
}

=item $cl->setSaturation $saturation 

Sets the $cl's saturation in the HSV model. Valid range [0 .. 1].

=cut 

sub setSaturation {
	my $this=shift @_;
	my $sat=shift @_;
	my ($h,$s,$v)=$this->asHSV;
	$this->setHSV($h,$sat,$v);
}

=item $cl->rotHue $degrees 

Rotates the $cl's hue in the HSV/L model. Valid range [-360 .. 360].

=cut 

sub rotHue {
	my $this=shift @_;
	my $rot=shift @_;
	my ($h,$s,$v)=$this->asHSV;
	$h+=$rot;
	$h%=360;
	$this->setHSV($h,$s,$v);
}

=item $cl->setHue $hue 

Sets the $cl's hue in the HSV/L model. Valid range [0 .. 360].

=cut 

sub setHue {
	my $this=shift @_;
	my $hue=shift @_;
	my ($h,$s,$v)=$this->asHSV;
	$this->setHSV($hue,$s,$v);
}

=item $cl->addBrightness $brightness 

Adds to the $cl's brightness in the HSV model. Valid range [-1 .. 1].

=cut 

sub addBrightness {
        my $this=shift @_;
        my $vol=shift @_;
        my ($h,$s,$v)=$this->asHSV;
        $this->setHSV($h,$s,$v+$vol);
}

=item $cl->setBrightness $brightness

Sets the $cl's brightness in the HSV model. Valid range [0 .. 1].

=cut

sub setBrightness {
	my $this=shift @_;
	my $v=shift @_;
	my ($h,$s)=$this->asHSV;
	$this->setHSV($h,$s,$v);
}

=item $cl->addLightness $lightness 

Adds to the $cl's lightness in the HSL model. Valid range [-1 .. 1].

=cut 

sub addLightness {
        my $this=shift @_;
        my $vol=shift @_;
        my ($h,$s,$v)=$this->asHSL;
        $this->setHSL($h,$s,$v+$vol);
}

=item $cl->setLightness $lightness

Sets the $cl's lightness in the HSL model. Valid range [0 .. 1].

=cut

sub setLightness {
	my $this=shift @_;
	my $l=shift @_;
	my ($h,$s)=$this->asHSL;
	$this->setHSL($h,$s,$l);
}

1;

__END__

=back

=head1 AUTHOR

Alfred Reibenschuh alfredreibenschuh@yahoo.com.

=head1 HISTORY

version 0.1_01 -- first public test release

=head1 BUGS

Some ... please report them.

=head1 TODO

more color spaces ?

=cut
