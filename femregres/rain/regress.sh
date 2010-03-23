#!/bin/sh

[ -f PASSED ] && rm -f PASSED

./splitext < ~/bin/CR

./regress.pl z.*

if [ $? -eq 0 ]; then
  echo "test passed"
  touch PASSED
fi
