#!/bin/sh
#
# computes running time and estimates time to complete from log file
#
# revision log :
#
# ......2017	ggu	started from scratch
# 05.03.2018	ggu	adapted to new write statement of shyfem
# 30.05.2018	ggu	more info to terminal (end date)
#
# still to do: command line options:
#
#	-running	show only running processes
#	-finished	show only finished processes
#
#---------------------------------------------------------

debug=0
if [ "$1" = "-debug" ]; then
  debug=1
  shift
fi

if [ $# -eq 0 ]; then
  echo "Usage: $0 [-debug] log-file(s)"
  exit 1
fi

server="nemunas.ku.lt"
server=$( hostname )

#---------------------------------------------------------

Convert_to_seconds()
{
  date -u -d "$1" +"%s"
}

Convert_to_iso()
{
  date -u -d "$1" +"%Y-%m-%d::%H:%M:%S"
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
    echo "error"
  fi
}

Get_last_line()
{
  last_line=$( tail -5 $file | head -1 )
  final_date=$( CheckLogfile "$last_line" )

  if [ "$final_date" = "error" ]; then		#still running
    finished="NO"
    perc=$( echo $last_line | sed -E 's/.*[0-9]\s+//' )	#look for % sign
    if [ "$perc" = "%" ]; then
      iter1=$( echo $last_line | sed -E 's/\s*\/.*\%//' | sed -E 's/.*\s//' )
      iter2=$( echo $last_line | sed -E 's/.*\/\s*//' | sed -E 's/\s.*//' )
      #echo "iter 1/2: $iter1 $iter2"
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
    end=$( Convert_to_iso "$final_date" )
    echo "$file:  total: $total  simulation has finished at $end UTC"
  fi
}

Estimate() 
{

file=$1

start_line=$( head -3 $file | tail -1 )
start_date=$( CheckLogfile "$start_line" )
if [ "$start_date" = "error" ]; then
  echo "$file : no log file..."
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

echo "$file:  total: $total  done: $done  todo: $todo"

}

#---------------------------------------------------------

for file in $*
do
  #echo "estimating $file"
  Estimate $file
done

#---------------------------------------------------------

