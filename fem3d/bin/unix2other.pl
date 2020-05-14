#!/usr/bin/perl -spi.bak
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# transforms unix fortran into other
#
# dos:		cdos#
# lahey:	clahey#


  s/^cdos\#// if $dos;
  s/^clahey\#// if $lahey;

  s/^\#dos\#// if $dos;
  s/^\#lahey\#// if $lahey;

