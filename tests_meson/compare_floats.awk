join  $1 $2 | awk '
       function abs(v) {return v < 0 ? -v : v} 
       function max(v1,v2) {return v1<v2 ? v2 : v1 } 

       BEGIN{
           max_err=0; 
           sum_err2=0;
           count=1;
           }
 
       ($1!=0){
           err=abs($2-$3); 
           max_err=max(max_err,err); 
           sum_err2=sum_err2+err*err ; 
           #printf("%e , %e, %e \n",err,max_err,sum_err2)
           }

       ($1==0 && NR!=1){
           printf("iter: %d diff error norms are L2: %e Linfty: %e\n",count,sqrt(sum_err2),max_err)
           max_err=0; 
           sum_err2=0;
           count+=1;
           }
       END{
           printf("iter: %d diff error norms are L2: %e Linfty: %e\n",count,sqrt(sum_err2),max_err)
           } 
     '
