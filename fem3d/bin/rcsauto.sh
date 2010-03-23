#!/bin/sh

for file
do
  echo $file
  ci -l -f $file << EOI
SHYFEM as of 26.06.1997 (with a lof of structural changes)
.
EOI
done

#original subroutines for FEM model as on 28.05.97
