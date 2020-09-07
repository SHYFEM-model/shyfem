#!/bin/bash
#
# examples and usage of advanced features of bash
#
#-----------------------------------------------------------

Loop()
{
  echo "----------------Loop over range---------------------"
  for id in {1..5};
  do
    echo "$id"
  done
}

Array()
{
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

  echo "--------------------Pop one item--------------------"
  allThreads=("${allThreads[@]:1}")
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
  Loop
  Array
  #Associative
}

Test

#-----------------------------------------------------------

