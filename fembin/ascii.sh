#!/bin/sh

grep --color='auto' -P -n "[\x80-\xFF]" $*

