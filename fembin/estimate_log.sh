#!/bin/sh
#
# computes running time and estimates time to complete from log file
#
#---------------------------------------------------------

if [ $# -eq 0 ]; then
  echo "Usage: $0 log-file"
  exit 1
fi

debug=1
debug=0

server="nemunas.ku.lt"

#---------------------------------------------------------

Convert_to_seconds()
{
  date -u -d "$1" +"%s"
}

Convertsecs() 
{
    h=$(($1/3600))
    m=$((($1/60)%60))
    s=$(($1%60))
    printf "%02d:%02d:%02d" $h $m $s
}

CheckLogfile()
{
  check=$( echo $1 | sed -e s'/.* - //' )
  #>&2 echo "check: $check"
  if [ "$check" = $server ]; then
    date=$( echo $1 | sed -e s'/ - .*//' )
    echo "$date"
  else
    #>&2 echo "cannot parse date line: $1"
    #echo "error"
    #try anyway...
    date=$( echo $1 | sed -e s'/ - .*//' )
    echo "$date"
  fi
}

Get_last_line()
{
  last_line=$( tail -5 $file | head -1 )
  final_date=$( CheckLogfile "$last_line" )

  if [ "$final_date" = "error" ]; then		#still running
    finished="NO"
    last_line=$( tail -1 $file )
    perc=$( echo $last_line | cut -d " " -f8 )	#look for % sign
    if [ "$perc" = "%" ]; then
      iter1=$( echo $last_line | cut -d " " -f4 )
      iter2=$( echo $last_line | cut -d " " -f6 )
    else
      finished="YES"
      echo "$file:  error determining iterations... run again"
    fi
  else
    finished="YES"
    start_secs=$( Convert_to_seconds "$start_date" )
    final_secs=$( Convert_to_seconds "$final_date" )
    total_secs=$(( $final_secs - $start_secs ))
    total=$( Convertsecs $total_secs )
    echo "$file  total: $total  simulation has finished..."
  fi
}

Estimate() 
{

file=$1

start_line=$( head -3 $file | tail -1 )
start_date=$( CheckLogfile "$start_line" )
if [ "$start_date" = "error" ]; then
  echo "$file : no log file or error reading date..."
  return
fi
end_date=$( date )

Get_last_line
[ "$finished" = "YES" ] && return

start_secs=$( Convert_to_seconds "$start_date" )
end_secs=$( Convert_to_seconds "$end_date" )
done_secs=$(( $end_secs - $start_secs ))
est_secs=$(( $iter2 * $done_secs / $iter1 ))
todo_secs=$(( $est_secs - $done_secs ))

if [ $debug -ne 0 ]; then
  echo "start: $start_date"
  echo "end:   $end_date"
  echo "last:  $last"
  echo "iterations: $iter1 - $iter2"

  echo "start secs: $start_secs"
  echo "end secs:   $end_secs"
  echo "done secs:  $done_secs"
  echo "est secs:   $est_secs"
  echo "todo secs:  $todo_secs"
fi

done=$( Convertsecs $done_secs )
todo=$( Convertsecs $todo_secs )
total=$( Convertsecs $est_secs )

echo "$file  total: $total  done: $done  todo: $todo"

}

#---------------------------------------------------------

for file in $*
do
  #echo "estimating $file"
  Estimate $file
done

#---------------------------------------------------------

