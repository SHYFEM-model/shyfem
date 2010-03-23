#!/usr/bin/perl
#
# this installs all modules in @mods or adjourns them
# it can be run anytime (in the worst case it does nothing)
# it must be run as root
# the CPAN module must be installed first by hand (see INSTALL)

use CPAN;

@mods = (
		 "MIME::Base64"
		,"Mail::Sender"
	);

#$mod = "Mail::Sender";
#install $mod;

foreach my $mod (@mods) {
    install $mod;
}

