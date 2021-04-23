#!/bin/bash
# takes the following arguments:
# 1. the number of mpi procs to be used
# 2. the executable to run
# 3.the input file 
# 4. the reference output file
# 5. the name of the generated output file
# 6. the test directory

_num_mpi_procs=$(($1))
_shyfem_exe=$2
_shyfem_input=$3
_ref_output_base_name=$4
_gen_output_base_name=$5
_test_dir=$6
ln -s ${_test_dir}mesh.${_num_mpi_procs}.bas mesh.${_num_mpi_procs}.bas > /dev/null 2>&1
ln -s ${_test_dir}wind.txt wind.txt > /dev/null 2>&1
if [ $_num_mpi_procs -gt 1 ]
then
  _shy_flag='-mpi'
else
  _shy_flag=''
fi

rm ${_gen_output_base_name}_* > /dev/null 2>&1
echo "running tests in dir $(pwd) : mpirun -np $_num_mpi_procs $_shyfem_exe $_shy_flag ${_test_dir}/$_shyfem_input > ${_shyfem_input}.out" 
mpirun -np $_num_mpi_procs $_shyfem_exe $_shy_flag ${_test_dir}/$_shyfem_input > ${_shyfem_input}.out 

# test if run was successfull 
if [[ $? -ne 0 ]]
then
    echo "running shyfem $_shyfem_input with $_num_mpi_procs mpi processes failed"
    exit 1
fi
#for n in {1..$_num_mpi_procs}

for (( n=1; n<=$_num_mpi_procs; n++ )); do
  # compare reference and new generated output
  _gen_output_full_name=${_gen_output_base_name}_$(printf "%04d" $n)
  _ref_output_full_name=${_test_dir}/${_ref_output_base_name}_$(printf "%04d" $n)
  diff $_ref_output_full_name $_gen_output_full_name > /dev/null 2>&1

  # test if diff was successfull
  if [[ $? -ne 0 ]]
  then
    echo "ERROR test failed running shyfem ${_test_dir}/$_shyfem_input with $_num_mpi_procs mpi processes : $_gen_output_full_name differs from $_ref_output_full_name"
    exit 1
  fi
done
rm  *.inf > /dev/null 2>&1
if [ $_num_mpi_procs -gt 1 ] ; then
  rm mpi_debug_* > /dev/null 2>&1
fi
exit 0


