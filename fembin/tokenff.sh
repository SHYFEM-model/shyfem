#!/bin/sh

nocom=NO

while [ -n "$1" ]
do
  case $1 in
	-nocom)		nocom=YES;;
	*)		break;;
  esac
  shift
done


if [ $nocom = YES ]; then
  token $* | grep -vi ': c' | grep -v '\s*!'
else
  token $*
fi


