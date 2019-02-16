#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

for file
do
  echo $file
  ci -l -f $file << EOI
SHYFEM as of 26.06.1997 (with a lof of structural changes)
.
EOI
done

#original subroutines for FEM model as on 28.05.97
