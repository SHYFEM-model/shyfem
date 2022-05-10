#!/bin/bash

if [ ! -f ggg ]; then
  echo "no file ggg... create first with make fem &> ggg"
  exit 1
fi

window="-B 4"

grep "Warning:" ggg | \
	grep -v "character" | \
	grep -v 'Nonstandard type declaration REAL' | \
	grep -v "Statement function at" | \
	grep -v "Computed GOTO a"

echo "============ specific ============"
grep $window "is used before it is typed at" ggg
grep $window "Old-style initialization at" ggg
grep $window "Nonstandard type declaration INTEGER" ggg
grep $window "to INTEGER" ggg
grep $window "to INTEGER" ggg | wc

