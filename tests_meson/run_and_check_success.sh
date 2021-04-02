#!/bin/bash
# takes the following arguments:
# 1. the number of mpi procs to be used
# 2. the executable to run
# 3.the input file 
# 4. the reference output file
# 5. the name of the generated output file
# 6. the test directory
# 7. the threshold value of the L2 and Linfty error-norm between ref and generated file

_num_mpi_procs=$(($1))
_shyfem_exe=$2
_shyfem_input=$3
_ref_output_base_name=$4
_gen_output_base_name=$5
_test_dir=$6
_threshold=$7
_petscrc_file=$8

if [ ! -f mesh.${_num_mpi_procs}.bas ]; then
  echo "ERROR : bas file mesh.${_num_mpi_procs}.bas does not exist"
  exit 1
fi

#ln -s ${_test_dir}mesh.${_num_mpi_procs}.bas mesh.${_num_mpi_procs}.bas > /dev/null 2>&1
if [ ! -f wind.txt ]; then
   ln -s ${_test_dir}/wind.txt wind.txt > /dev/null 2>&1
fi
if [ ! -f "AmgX.info" ]; then
   cp ${_test_dir}/"AmgX.info" .  > /dev/null 2>&1
fi

if [ ! -f $_shyfem_input ]; then
   echo "cp ${_test_dir}/$_shyfem_input "
   cp ${_test_dir}/$_shyfem_input .  
fi

if [ $_petscrc_file != 'none' ] ; then
  if [ ! -f $_petscrc_file ]; then
     cp ${_test_dir}/$_petscrc_file .  > /dev/null 2>&1
  fi
  _petsc_config="-zeta_petsc_config $_petscrc_file"
  _shyfem_output=$( basename ${_petscrc_file} .petscrc ).out
else
  _petsc_config=' ' 
  _shyfem_output=$( basename ${_shyfem_input} .str ).out
fi

if [ $_num_mpi_procs -gt 1 ]
then
  _shy_flag='-mpi'
else
  _shy_flag=''
fi

rm ${_shyfem_output} > /dev/null 2>&1
rm ${_gen_output_base_name}_* > /dev/null 2>&1
echo "running tests in dir $(pwd) : mpirun -np $_num_mpi_procs $_shyfem_exe $_shy_flag $_petsc_config $_shyfem_input > ${_shyfem_output}" 
mpirun -np $_num_mpi_procs $_shyfem_exe $_shy_flag  $_petsc_config $_shyfem_input > ${_shyfem_output} 
returned_val=$?
normal_exit_string=$( grep -q "shyfem program exiting normally" ${_shyfem_output} ; echo $?)
# test if run was successfull 
if [[ $returned_val -ne 0 ]] || [ $normal_exit_string -ne 0 ]
then
    echo "running shyfem $_shyfem_input with $_num_mpi_procs mpi processes failed"
    exit 1
fi
if [ _threshold -ne 'none' ] ; then
  _threshold=$(($7))
  #for n in {1..$_num_mpi_procs}
  for (( n=1; n<=$_num_mpi_procs; n++ )); do
    # compare reference and new generated output
    _gen_output_full_name=${_gen_output_base_name}_$(printf "%04d" $n)
    _ref_output_full_name=${_test_dir}/${_ref_output_base_name}_$(printf "%04d" $n)

    RESULTS=$(join $_ref_output_full_name $_gen_output_full_name | awk -v threshold=$_threshold '
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
                       #printf("  --- %e,%e, %e , %e, %e \n",$2,$3,err,max_err,sum_err2)
                       }

                   ($1==0 && NR!=1){
                       if(sqrt(sum_err2)>=threshold || max_err>=threshold)
                         printf("  >>> iter: %d diff error norms are L2: %e Linfty: %e\n",count,sqrt(sum_err2),max_err)
                       max_err=0; 
                       sum_err2=0;
                       count+=1;
                       }
                   END{
                       if(sqrt(sum_err2)>=threshold || max_err>=threshold)
                          printf("  >>> iter: %d diff error norms are L2: %e Linfty: %e\n",count,sqrt(sum_err2),max_err)
                       } 
                 ' ) 
                 #| awk '($8>1E-08 || $10>1E-08){printf("%e %e \n",$8,$10)}' )

    if [[ $(echo -n $RESULTS | wc -c)>0 ]] ; then
      echo " " 
      echo "ERROR some L2 or Linfty error norm between ref and generated files are equals or greater than threshold $_threshold : "
      echo "$RESULTS"
      echo "floating point comparison test FAILED running shyfem $_shyfem_input with $_num_mpi_procs mpi processes : $_gen_output_full_name differs from $_ref_output_full_name"
      echo " " 
      exit 1
    fi

  done
fi
rm  *.inf > /dev/null 2>&1
if [ $_num_mpi_procs -gt 1 ] ; then
  rm mpi_debug_* > /dev/null 2>&1
fi
exit 0


