#!/usr/bin/perl -spi.bak
#
# $Id: unix2other.pl,v 1.2 1999/02/12 15:34:35 georg Exp $
#
# transforms unix fortran into other
#
# dos:		cdos#
# lahey:	clahey#


  s/^cdos\#// if $dos;
  s/^clahey\#// if $lahey;

  s/^\#dos\#// if $dos;
  s/^\#lahey\#// if $lahey;

