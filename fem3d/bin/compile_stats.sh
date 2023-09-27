
wc LIST*

old=""

for list in LIST[1-9]
do
  if [ -n "$old" ]; then
    echo "checking $old $list"
    cat $old $list | sort | uniq -u
  fi
  old=$list
done

