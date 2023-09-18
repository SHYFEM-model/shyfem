#!/bin/bash

for FILE in `ls *.f*`; 
   do 
      VAR=`md5sum $FILE | cut -d\  -f1`
      #echo $VAR
      VAR2=`md5sum /work/asc/gm30419/shympi/test/merge/11_09/src/$FILE | cut -d\  -f1`
      #echo $VAR $VAR2
      if [ "$VAR" != "$VAR2" ] ; then
         echo "error $FILE $VAR $VAR2 "
      else
         echo $FILE $VAR ok
      fi
   done
