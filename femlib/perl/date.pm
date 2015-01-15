#!/usr/bin/perl -w
#
# date utilities
#
# example of usage:
#
# use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
#
# use date;
# 
# my $date = new date;
#
# $date->init_year($year0);		#set reference year
# $it = $date->convert_to_it($year,$month,$day,$hour,$min,$sec);
# ($year,$month,$day,$hour,$min,$sec) = $date->convert_from_it($it);
#
# to convert date online please use the program convert_date.pl
#
##############################################################
#
# 10.03.2010	ggu	bug in adjust_mod() -> account for negative values
# 28.09.2010	ggu	new test routine test_femdate()
# 21.10.2014	ggu	absolute time routines introduced
# 21.10.2014	ggu	refer to 1.1.1 (date2days and days2date)
# 05.11.2014	ggu	new routine unformat_time_date()
#
# version 2.2
#
##############################################################

use strict;

package date;

##############################################################

sub new
{
    my $self;

    $self =	{
	    		 days0		=>	undef
			,secs0		=>	undef
			,verbose	=>	0
			,md		=> 	[]
			,td		=> 	[]
			,daynames	=>	{}
			,monthnames	=>	{}
			,it_daynames	=>	{}
			,it_monthnames	=>	{}
		};

    bless $self;
    $self->init_arrays();
    $self->init_it(19700101,0);		# 1/1/1970 is it=0
    return $self;
}

###############################################################################

sub init_arrays
{
	my ($self) = @_;

	my @m = qw/0 31 28 31 30 31 30 31 31 30 31 30 31/;
	my @t = ();
	my @skew = qw/0 0 3 3 6 1 4 6 2 5 0 3 5/;

	my @days = qw/mon tue wed thu fri sat sun/;
	my @months = qw/jan feb mar apr may jun jul aug sep oct nov dec/;
	my @it_days = qw/lun mar mer gio ven sab dom/;
	my @it_months = qw/gen feb mar apr mag giu lug ago set ott nov dic/;

	my $tot = 0;
	for(my $i=0;$i<=12;$i++) {
	  $tot += $m[$i];
	  $t[$i] = $tot;
	}

	$self->{md} = \@m;
	$self->{td} = \@t;
	$self->{skew} = \@skew;
	$self->{daynames} = $self->make_names_to_numbers(\@days);
	$self->{monthnames} = $self->make_names_to_numbers(\@months);
	$self->{it_daynames} = $self->make_names_to_numbers(\@it_days);
	$self->{it_monthnames} = $self->make_names_to_numbers(\@it_months);
}

sub make_names_to_numbers
{
	my ($self,$list) = @_;
	my %hash = ();

	my $i = 0;
	foreach my $key (@$list) {
	  $hash{$key} = ++$i;
	}

	return \%hash;
}

sub bises
{
	my ($self,$y) = @_;

	if( ($y%4 == 0 and $y%100 != 0) or $y%400 == 0 ) {
	  return 1;
	} else {
	  return 0;
	}
}

sub total_days_in_year
{
	my ($self,$year) = @_;

	if( $self->bises($year) ) {
	  return 366;
	} else {	
	  return 365;
	}
}

sub days_in_month
{
	my ($self,$year,$month) = @_;

	my $days = $self->{md}->[$month];
	$days++ if $month == 2 and $self->bises($year);

	return $days;
}

sub total_days_in_month
{
	my ($self,$year,$month) = @_;

	my $days = $self->{td}->[$month];
	$days++ if $month >= 2 and $self->bises($year);

	return $days;
}

sub from_daynames 
{
	my ($self,$dayname) = @_;

	my $day;
	$dayname = lc($dayname);

	$day = $self->{daynames}->{$dayname};
	return $day if $day;

	$day = $self->{it_daynames}->{$dayname};
	return $day if $day;

	return "";
}

sub from_monthnames
{
	my ($self,$monthname) = @_;

	my $month;
	$monthname = lc($monthname);

	$month = $self->{monthnames}->{$monthname};
	return $month if $month;

	$month = $self->{it_monthnames}->{$monthname};
	return $month if $month;

	return "";
}

#--------------------------------------

sub convert_args		# internal routine
{
	my $self = shift;

	my ($year,$month,$day);
	my ($hour,$min,$sec);

	my @args = @_;

	$args[0] = 0 unless $args[0];
	$args[1] = 0 unless $args[1];

	my $date = int($args[0]+0.5);	#convert to numeric
	my $time = int($args[1]+0.5);

	if( $date > 10000 ) {		# must convert
	  ($year,$month,$day) = $self->unpack_date($date);
	  ($hour,$min,$sec)   = $self->unpack_time($time);
	} else {			# already converted
	  ($year,$month,$day,$hour,$min,$sec) = @_;
	  $month = 1 unless $month;
	  $day = 1 unless $day;
	  $hour = 0 unless $hour;
	  $min = 0 unless $min;
	  $sec = 0 unless $sec;
	}

	return ($year,$month,$day,$hour,$min,$sec);
}

sub init_year		#helper function
{
	my ($self,$year) = @_;

	$self->init_it($year,1,1,0,0,0);
}

sub init_date		#helper function
{
	my ($self,$date) = @_;

	my ($year,$month,$day);
	my ($hour,$min,$sec) = (0,0,0);

	if( $date =~ /^\'(.*)\'$/  ) {		#string
	  $date = $1;
	  ($year,$month,$day,$hour,$min,$sec)=$self->unformat_time_date($date);
	} elsif( $date =~ /:/  ) {		#string
	  ($year,$month,$day,$hour,$min,$sec)=$self->unformat_time_date($date);
	} elsif( $date < 10000 ) {
	  $year = $date; $month = 1; $day = 1;
	} else {
	  ($year,$month,$day) = $self->unpack_date($date);
	}

	$self->init_it($year,$month,$day,$hour,$min,$sec);
}

sub init_it
{
	my $self = shift;

	my ($year,$month,$day,$hour,$min,$sec) = $self->convert_args(@_);

	my $days0 = $self->date2days($year,$month,$day);
	my $secs0 = $self->time2secs($hour,$min,$sec);

	#print "init_it: $_[0]  $_[1]\n";
	#print "init_it: $year,$month,$day,$hour,$min,$sec\n";
	#print "init_it: $days0 $secs0\n";

	$self->{days0} = $days0;
	$self->{secs0} = $secs0;
}

#--------------------------------

sub convert_to_abs
{
	my $self = shift;

	my ($year,$month,$day,$hour,$min,$sec) = $self->convert_args(@_);

	my $days = $self->date2days($year,$month,$day);
	my $secs = $self->time2secs($hour,$min,$sec);

	return 86400.*$days + $secs;
}

sub convert_from_abs
{
	my ($self,$atime) = @_;

        my $days = int($atime/86400.);
        my $secs = $atime - 86400.*$days;

        my ($year,$month,$day) = $self->days2date($days);
        my ($hour,$min,$sec) = $self->secs2time($secs);

        my $date = $self->pack_date($year,$month,$day);
        my $time = $self->pack_time($hour,$min,$sec);

	return ($date,$time);
}

sub format_abs
{
	my ($self,$atime) = @_;

	my ($date,$time) = $self->convert_from_abs($atime);
	return $self->format_time_date($date,$time);
}

#--------------------------------

sub convert_to_it
{
	my $self = shift;

	my ($year,$month,$day,$hour,$min,$sec) = $self->convert_args(@_);

	my $days = $self->date2days($year,$month,$day);
	my $secs = $self->time2secs($hour,$min,$sec);

	my $days0 = $self->{days0};
	my $secs0 = $self->{secs0};

	return ($days-$days0)*86400 + ($secs-$secs0);
}

sub convert_from_it
{
	my ($self,$it) = @_;

	my ($year,$month,$day) = $self->days2date($self->{days0});
	my ($hour,$min,$sec)   = $self->secs2time($self->{secs0});

	$sec += $it;
	
	($sec,$min) = $self->adjust_mod($sec,$min,60);
	($min,$hour) = $self->adjust_mod($min,$hour,60);
	($hour,$day) = $self->adjust_mod($hour,$day,24);

	#print STDERR "**** $sec,$min    $min,$hour    $hour,$day\n";
	# we go back to first of year

	my $jday = $self->date2julian($year,$month,1);
	$day += $jday - 1;
	$month = 1;

	# get positive days

	while( $day <= 0 ) {
	  $year--;
	  $day += $self->total_days_in_year($year);
	}

	# go to year
	
	while( $day > $self->total_days_in_year($year) ) {
	  $day -= $self->total_days_in_year($year);
	  $year++;
	}

	# adjust days finally
	
	$jday = $day;
	($month,$day) = $self->julian2date($jday,$year);

	return ($year,$month,$day,$hour,$min,$sec);
}

sub adjust_mod			# internal routine
{
	my ($self,$low,$high,$mod) = @_;

	return ($low,$high) if 0 <= $low and $low < $mod;

	my $ia = int($low/$mod);

	$low -=  $ia*$mod;
	$high += $ia;

	if( $low < 0 ) {		#bug fix 10.03.2010
	  $low += $mod;
	  $high -= 1;
	}

	die "Cannot convert time: $low $high\n" if $low < 0;

	return ($low,$high);
}

#--------------------------------------

# next algorithm computes weekday from date 
# 0 = sunday, 1 = monday, etc.. 6 = saturday
# for algorithm see http://klausler.com/new-dayofweek.html

sub weekday
{
	my ($self,$year,$month,$day) = @_;

	die "Cannot use year > 2099 in weekday\n" if $year > 2099;

	my $sum = $year;
	my $skew = $self->{skew};

	$sum -= 1900;
	$sum += int($sum/4);			# years skew value
	$sum-- if $self->bises($year) and $month <= 2;
	$sum += $day;
	$sum += $skew->[$month];		# months skew value

	return $sum % 7;
}

#--------------------------------------

sub date2julian
{
	my ($self,$year,$month,$day) = @_;

	return $day + $self->total_days_in_month($year,$month-1);
}

sub julian2date
{
	my ($self,$jd,$year) = @_;

	my ($month,$day);

	$month = int($jd/30);
	$month-- if $jd <= $self->total_days_in_month($year,$month);
	$day = $jd - $self->total_days_in_month($year,$month);
	$month++;

	return ($month,$day);
}

#--------------------------------------

sub date2days		#days from 1/1/1
{
	my ($self,$year,$month,$day) = @_;

	if( $month < 1 or $month > 12 ) {
	  die "** error: month must be in [1-12]: $month\n";
	}
	if( $day < 1 or $day > 31 ) {
	  die "** error: day must be in [1-31]: $day\n";
	}

	my $y = $year - 1;
	my $days = 365*$y + int($y/4) - int($y/100) + int($y/400);

	my $jd = $self->date2julian($year,$month,$day);
	#return $days + $jd;	#old - one day off
	return $days + $jd - 1;
}

sub days2date
{
	my ($self,$days) = @_;

	#my $days1 = $days - 1;	#old - one day off
	my $days1 = $days;

	my $y = 2 + int( (400*$days)/(365*400+97) );	#rough first estimate

	my $aux = $days1 + 1;
	while( $aux > $days1 ) {
	  $y--;
	  $aux = 365*$y + int($y/4) - int($y/100) + int($y/400);
	}

	my $year = $y + 1;
	$aux = $days1 - $aux + 1;

	my ($month,$day) = $self->julian2date($aux,$year);

	return ($year,$month,$day);
}

sub time2secs
{
	my ($self,$hour,$min,$sec) = @_;

	return $sec + 60*$min + 3600*$hour;
}

sub secs2time
{
	my ($self,$secs) = @_;

	my ($hour,$min,$sec);

	$min = int($secs/60);
	$sec = $secs - $min*60;

	$hour = int($min/60);
	$min = $min - $hour*60;

	return ($hour,$min,$sec);
}

#--------------------------------------

sub pack_date
{
	my ($self,$year,$month,$day) = @_;

	return $day + 100*$month + 10000*$year;
}

sub unpack_date
{
	my ($self,$date) = @_;

	my ($year,$month,$day);

	$month = int($date/100);
	$day = $date - $month*100;

	$year = int($month/100);
	$month = $month - $year*100;

	return ($year,$month,$day);
}

sub pack_time
{
	my ($self,$hour,$min,$sec) = @_;

	return $sec + 100*$min + 10000*$hour;
}

sub unpack_time
{
	my ($self,$time) = @_;

	my ($hour,$min,$sec);

	$min = int($time/100);
	$sec = $time - $min*100;

	$hour = int($min/100);
	$min = $min - $hour*100;

	return ($hour,$min,$sec);
}

sub format_time_date
{
        my $self = shift;

        my ($year,$month,$day,$hour,$min,$sec) = $self->convert_args(@_);

	$month = $self->add_zero2($month);
	$day = $self->add_zero2($day);
	$hour = $self->add_zero2($hour);
	$min = $self->add_zero2($min);
	$sec = $self->add_zero2($sec);

	my $line = "$year-$month-${day}::$hour:$min:$sec";

	return $line;
}

sub unformat_time_date
{
        my ($self,$line) = @_;

	my ($date,$time) = split("::",$line);
	my ($year,$month,$day) = split("-",$date);
	my ($hour,$min,$sec) = split(":",$time);

	$year = 0 unless $year;
	$month = 0 unless $month;
	$day = 0 unless $day;
	$hour = 0 unless $hour;
	$min = 0 unless $min;
	$sec = 0 unless $sec;

        return ($year,$month,$day,$hour,$min,$sec);
}

sub add_zero2
{
        my ($self,$val) = @_;

	$val = "0$val" if( $val < 10 and not $val =~ /^0/ );
	$val = "00"    if( $val == 0 );

	return $val;
}

#----------------------------------------------------------------

sub parse_email_date
{
        my ($self,$line) = @_;

	my $error = 0;
	my ($wd,$day,$month,$year,$time);

# p1	Fri, 22 Jun 2001 10:32:34 +0200
# p2	22 Jun 2001 10:32:34 +0200
# p3	Mon Dec 1 17:33:37 CST 1997
# p4	Mon Jun 18 13:55:01 2001

	my $p1 = '(\w{3}),\s+(\d+)\s+(\w{3})\s+(\d+)\s+([0-9:\.]+)(.*)';
	my $p2 = '(\d+)\s+(\w{3})\s+(\d+)\s+([0-9:\.]+)(.*)';
	my $p3 = '(\w{3})\s+(\w{3})\s+(\d+)\s+([0-9:\.]+)\s+(\w{3})\s+(\d+)';
	my $p4 = '(\w{3})\s+(\w{3})\s+(\d+)\s+([0-9:\.]+)\s+(\d+)';

	return 0 unless $line;

	if( $line =~ /^$p1$/ ) {
	  ($wd,$day,$month,$year,$time) = ($1,$2,$3,$4,$5);
	} elsif( $line =~ /^$p2$/ ) {
	  ($wd,$day,$month,$year,$time) = ("Sun",$1,$2,$3,$4);
	} elsif( $line =~ /^$p3$/ ) {
	  ($wd,$day,$month,$year,$time) = ($1,$3,$2,$6,$4);
	} elsif( $line =~ /^$p4$/ ) {
	  ($wd,$day,$month,$year,$time) = ($1,$3,$2,$5,$4);
	} else {
	  $error = 1;
	}

	if( $error ) {
	  print STDERR "*** Cannot parse date: |$line|\n";
	} else {
	  ($day,$month,$year,$error) = $self->check_date($wd,$day,$month,$year);
	  $error = $self->check_time($time) unless $error;
	  if( $error ) {
	    print STDERR "*** Error parsing date: |$line|\n";
	  }
	}

	if( $error ) {
	  exit 1;
	} else {
	  my $itdate = $self->convert_to_it($year,$month,$day,0,0,0);
	  #print STDERR "Date: $line -> $year-$month-$day ($itdate)\n";
	  return $itdate;
	}
}

sub check_time
{
        my ($self,$time) = @_;

	my $error = 0;
	my ($hour,$min,$sec);

	if( $time =~ /\d{1,2}:\d{1,2}:\d{1,2}/ ) {
	  ($hour,$min,$sec) = ($1,$2,$3);
	  $error = 0;
	} elsif( $time =~ /\d{1,2}:\d{1,2}/ ) {
	  ($hour,$min,$sec) = ($1,$2,0);
	  $error = 0;
	} elsif( $time =~ /\d{1,2}\.\d{1,2}\.\d{1,2}/ ) {
	  ($hour,$min,$sec) = ($1,$2,$3);
	  $error = 0;
	} else {
	  $error = 1;
	}

	if( $error ) {
	  print STDERR "*** error parsing time: $time\n";
	} else {
	  $error = 1;
	  if( $hour < 0 or $hour >= 24 ) {
	    print STDERR "*** error in hour: $hour  $time\n";
	  } elsif( $min < 0 or $min >= 60 ) {
	    print STDERR "*** error in minute: $min  $time\n";
	  } elsif( $sec < 0 or $sec >= 60 ) {
	    print STDERR "*** error in second: $sec  $time\n";
	  } else {
	    $error = 0;
	  }
	}

	return $error;
}

sub check_date
{
        my ($self,$wd,$day,$month,$year) = @_;

	my $wd1 = $self->from_daynames($wd);
	my $month1 = $self->from_monthnames($month);

	my $error = 1;

	if( $wd1 < 1 ) {
	  print STDERR "*** error parsing weekday: $wd\n";
	} elsif( $day < 1 or $day > 31 ) {
	  print STDERR "*** error parsing day: $day\n";
	} elsif( $month1 < 1 or $month1 > 12 ) {
	  print STDERR "*** error parsing month: $month\n";
	} elsif( $year < 1980 or $year > 2020 ) {
	  if( $year > 69 and $year <= 100 ) {
	    $year += 1900;
	    $error = 0;
	    #print STDERR "* Changed year... $year\n";
	  } elsif( $year <= 20 ) {
	    $year += 2000;
	    $error = 0;
	    #print STDERR "* Changed year... $year\n";
	  } else {
	    print STDERR "*** error parsing year: $year\n";
	  }
	  if( $year < 1980 or $year > 2020 ) {
	    print STDERR "*** Cannot adjust year: $year\n";
	    $error = 1;
	  }
	} else {
	  $error = 0;
	}

	$month = $self->add_zero2($month1);

	return ($day,$month,$year,$error);
}

##########################################################################

sub rand_num
{
	my ($self,$n) = @_;

	my $r = rand();

	my $i = 1 + int($r*$n);

	if( $i < 1 or $i > $n ) {
	  die "error in rand_num: $i $n\n";
	}
	#$i-- if $i > $n;

	return $i;
}

sub get_rand_date
{
	my ($self,$maxyear) = @_;

	my $year = $self->rand_num($maxyear);
	my $tdays = $self->total_days_in_year($year);
	my $jday = $self->rand_num($tdays);
	my ($month,$day) = $self->julian2date($jday,$year);

	return ($year,$month,$day,$jday);
}

sub test_equal
{
	my ($self,$i1,$i2) = @_;

	if( $i1 == $i2 ) {
	  return 0;
	} else {
	  return 1;
	}
}

sub test
{
	my ($self,$n) = @_;

	my $error = 0;

	$error = $self->test_julian($n);
	if( $error ) {
	  die "*** there were $error errors in test_julian\n";
	} else {
	  print STDERR "test_julian ok ($n)\n";
	}

	$error = $self->test_days($n);
	if( $error ) {
	  die "*** there were $error errors in test_days\n";
	} else {
	  print STDERR "test_days ok ($n)\n";
	}

	$error = $self->test_it($n);
	if( $error ) {
	  die "*** there were $error errors in test_it\n";
	} else {
	  print STDERR "test_it ok ($n)\n";
	}
}

#-----------------------------------------------------------

sub test_it
{
	my ($self,$n) = @_;

	my $error = 0;

	$self->init_it(20000101,0);

	$error += $self->test_random_it($n);

	$error += $self->test_single_it(86400);

	return $error;
}

sub test_random_it
{
	my ($self,$n) = @_;

	$n = 100 unless $n;
	my $error = 0;

	my $max = 10*365*86400;

	while( $n-- ) {
	  my $it = $self->rand_num($max) - $max/2;
	  $error += $self->test_single_it($it);
	}
}

sub test_single_it
{
	my ($self,$it) = @_;

	my ($year,$month,$day,$hour,$min,$sec) = $self->convert_from_it($it);
	my $it2 = $self->convert_to_it($year,$month,$day,$hour,$min,$sec);

	my $line = "";
	my $error = 0;
	if( $self->test_equal($it,$it2) ) {
	  $error = 1;
	  $line = "  ******";
	}
	if( $self->{verbose} ) {
	  print "$year $month $day $hour $min $sec  $it $it2   $line \n";
	}

	return $error;
}

#-----------------------------------------------------------

sub test_days
{
	my ($self,$n) = @_;

	my $error = 0;

	$error += $self->test_random_days($n);

	$error += $self->test_single_days(1);
	$error += $self->test_single_days(365);
	$error += $self->test_single_days(366);

	return $error;
}

sub test_random_days
{
	my ($self,$n) = @_;

	$n = 100 unless $n;
	my $error = 0;

	while( $n-- ) {
	  my $days = $self->rand_num(3000*365);
	  $error += $self->test_single_days($days);
	}

	return $error;
}

sub test_single_days
{
	my ($self,$days) = @_;

	my ($year,$month,$day) = $self->days2date($days);
	my $days2 = $self->date2days($year,$month,$day);

	my $line = "";
	my $error = 0;
	if( $self->test_equal($days,$days2) ) {
	  $error = 1;
	  $line = "  ******";
	}
	print "$year $month $day  $days $days2  $line \n" if $self->{verbose};

	return $error;
}

#-----------------------------------------------------------

sub test_julian
{
	my ($self,$n) = @_;

	my $error = 0;

	$error += $self->test_rand_julian($n);

	$error += $self->test_single_julian(1,1999,1,1);
	$error += $self->test_single_julian(31,1999,1,31);
	$error += $self->test_single_julian(365,1999,12,31);
	$error += $self->test_single_julian(366,2004,12,31);

	return $error;
}

sub test_single_julian
{
	my ($self,$jday,$year,$month,$day) = @_;

	my $jday2 = $self->date2julian($year,$month,$day);

	my $line = "";
	my $error = 0;
	if( $self->test_equal($jday,$jday2) ) {
	  $error = 1;
	  $line = "  ******";
	}
	print "$year $month $day  $jday $jday2  $line \n" if $self->{verbose};

	return $error;
}

sub test_rand_julian
{
	my ($self,$n) = @_;

	$n = 100 unless $n;
	my $error = 0;

	while( $n-- ) {
	  my ($year,$month,$day,$jday) = $self->get_rand_date(3000);
	  $error += $self->test_single_julian($jday,$year,$month,$day);
	}

	return $error;
}

sub test_weekday
{
	my $date = new date;

	while( 1 ) {
	  print "Please enter date (YYYYMMDD): ";
	  my $line = <>;
	  chomp($line);
	  last unless $line;
	  my ($year,$month,$day) = $date->unpack_date($line);
	  my $wd = $date->weekday($year,$month,$day);
	  print "($line)  $year $month $day    $wd\n";
	}
}

sub test_femdate
{
	my $date = new date;

	while( 1 ) {
	  print "Please enter date (YYYYMMDD): ";
	  my $line = <>;
	  chomp($line);
	  last unless $line;
	  my ($year,$month,$day) = $date->unpack_date($line);
	  $date->init_year($year);
	  my $it = $date->convert_to_it($year,$month,$day,0,0,0);
	  print "($line)  $year $month $day    $it\n";
	}
}

sub test_abs_time
{
	my $pdate = new date;

	my $niter = 100;
	my $date = 0;
	my $time = 0;
	my $line = "";

	my $atime = 86400;
	$line = $pdate->format_abs($atime);
	print "$line   $atime\n";

	$date = 1;
	$atime = $pdate->convert_to_abs($date,$time);
	$line = $pdate->format_abs($atime);
	print "$line   $atime\n";

	$date = 3000;
	$atime = $pdate->convert_to_abs($date,$time);
	$line = $pdate->format_abs($atime);
	print "$line   $atime\n";

	my $amax = $atime;
	my $dt = int($amax/$niter);

	for(my $i=1;$i<=$niter;$i++) {
	  $atime = $i*$dt;
	  ($date,$time) = $pdate->convert_from_abs($atime);
	  my $atime1 = $pdate->convert_to_abs($date,$time);
	  $line = $pdate->format_abs($atime);
	  print "  $atime  $atime1   $line  $date $time\n";
	  if( $atime != $atime1 ) {
	    my $adiff = $atime-$atime1;
	    die "*** $atime  $atime1   $adiff   $line\n";
	  }
	}
	
}


##########################################################################

sub test_this
{
	my $date = new date;

	print STDERR "running auto test...\n";
	$date->{verbose} = 0;
	$date->test(100);
}

#test_this() if $0 =~ /date.pm$/;
#test_weekday() if $0 =~ /date.pm$/;
test_abs_time() if $0 =~ /date.pm$/;

################################
1;
################################

