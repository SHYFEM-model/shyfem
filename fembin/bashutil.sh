#!/bin/bash
#
# examples and usage of advanced features of bash
#
#-----------------------------------------------------------

Associative()
{
  declare -A assArray1
  assArray1[fruit]=Mango
  assArray1[bird]=Cockatail
  assArray1[flower]=Rose
  assArray1[animal]=Tiger

  echo ${assArray1[bird]}
  echo ${assArray1[flower]}

  for key in "${!assArray1[@]}"; do echo $key; done
  echo "${!assArray1[@]}"

  for val in "${assArray1[@]}"; do echo $val; done
  echo "${assArray1[@]}"

  declare -A assArray2
  assArray2[monitor]=Dell
  assArray2[keyboard]=Samsung
  assArray2[graphic]=A4Tech

  echo "${assArray2[@]}"
  assArray2+=([Mouse]=Logitech)
  echo "${assArray2[@]}"

  unset assArray2[Monitor]
  echo ${assArray2[Monitor]}

  if [ ${assArray2[Monitor]+_} ]; then echo "Found"; else echo "Not found"; fi

  echo "${assArray1[@]}"
  unset assArray1
  echo "${assArray1[@]}"
}

log()
{ 
  # returns floor[log(x)] for x integer
  # call as "log 777" (base 2) or "log 10 777" (base 10)
  # http://phodd.net/gnu-bc/bcfaq.html#bashlog

  local x=$1 n=2 l=-1;
  if [ "$2" != "" ]; then 
    n=$x;x=$2
  fi

  while((x))
  do
    let l+=1 x/=n
  done

  echo $l;
}

#-----------------------------------------------------------

Test()
{
  Associative
}

Test

#-----------------------------------------------------------

