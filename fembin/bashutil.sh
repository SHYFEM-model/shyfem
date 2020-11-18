#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# examples and usage of advanced features of bash
#
#------------------------------------------------------------------------

func()
{
  return $1
}

If()
{
  echo "================================="
  echo "test if"
  echo "================================="

  if func 0; then echo "ok"; else echo "false"; fi
  if func 1; then echo "ok"; else echo "false"; fi
}

Loop()
{
  echo "================================="
  echo "test loop"
  echo "================================="

  echo "----------------Loop over range---------------------"
  for id in {1..5};
  do
    echo "$id"
  done
}

Array()
{
  echo "================================="
  echo "test array"
  echo "================================="

  MyArray=( HTML Javascript CSS JQuery Bootstrap )
 
  echo "----------Print 5 values individually---------------"
  echo ${MyArray[0]}
  echo ${MyArray[1]}
  echo ${MyArray[2]}
  echo ${MyArray[3]}
  echo ${MyArray[4]}
 
  echo "-----------------Print all values-------------------"
  echo ${MyArray[*]}
  echo "-----------------Print all values-------------------"
  echo ${MyArray[@]}

  echo "----------------Loop through values-----------------"
  for i in ${!MyArray[@]}; do
    echo "$i - ${MyArray[$i]}"
  done

  echo "--------------------Add to array--------------------"
  allThreads=(1 2 4 8 16 32 64 128)
  allRuntimes=()
  for t in ${allThreads[@]}; do
    runtime=$t
    allRuntimes+=( $runtime )
  done
  echo ${allThreads[*]}
  echo ${allRuntimes[*]}
  echo "array size: ${#allThreads[@]}"

  echo "------------------Shift first item------------------"
  allThreads=("${allThreads[@]:1}")
  echo ${allThreads[*]}
  echo "array size: ${#allThreads[@]}"

  echo "-------------------Pop last item--------------------"
  allThreads=("${allThreads[@]:0:${#allThreads[@]}-1}")
  echo ${allThreads[*]}
  echo "array size: ${#allThreads[@]}"

  echo "------------------Get items 2,3,4-------------------"
  allThreads=("${allThreads[@]:2:3}")
  echo ${allThreads[*]}
  echo "array size: ${#allThreads[@]}"

  # arr=(1 2 3) 	# Initialize array
  # ${arr[2]} 		# Retrieve third element
  # ${arr[@]} 		# Retrieve all elements
  # ${!arr[@]} 		# Retrieve array indices
  # ${#arr[@]} 		# Calculate array size
  # arr[0]=3 		# Overwrite 1st element
  # arr+=(4) 		# Append value(s)
  # str=$(ls) 		# Save ls output as a string
  # arr=( $(ls) ) 	# Save ls output as an array of files
  # ${arr[@]:s:n} 	# Retrieve n elements starting at index s
}

Associative()
{
  echo "================================="
  echo "test associative array"
  echo "================================="

  declare -A assArray1

  assArray1[fruit]=Mango
  assArray1[bird]=Cockatail
  assArray1[flower]=Rose
  assArray1[animal]=Tiger

  echo "------------------show keys of array-------------------"
  echo "${!assArray1[@]}"
  for key in "${!assArray1[@]}"
  do 
    echo $key
  done

  echo "------------------show vals of array-------------------"
  echo "${assArray1[@]}"
  for val in "${assArray1[@]}"
  do 
    echo $val
  done

  echo "------------------show keys and vals of array-------------------"
  for key in "${!assArray1[@]}"
  do 
    echo "$key: ${assArray1[$key]}"
  done

  echo "------------------access elements-------------------"
  echo "bird: ${assArray1[bird]}"
  echo "flower: ${assArray1[flower]}"

  echo "------------------add new element-------------------"
  echo "${assArray1[@]}"
  assArray1+=([vegatable]=Potato)
  echo "${assArray1[@]}"

  echo "------------------delete element-------------------"
  echo "${assArray1[@]}"
  unset assArray1[fruit]
  echo "${assArray1[@]}"
  echo "fruit: ${assArray1[fruit]}"

  echo "------------------look for element-------------------"
  echo -n "look for flower: "
  if [ ${assArray1[flower]+_} ]; then echo "Found"; else echo "Not found"; fi
  echo -n "look for fruit: "
  if [ ${assArray1[fruit]+_} ]; then echo "Found"; else echo "Not found"; fi

  echo "------------------delete array-------------------"
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
  Loop
  Array
  Associative
  If
}

Test

#-----------------------------------------------------------

