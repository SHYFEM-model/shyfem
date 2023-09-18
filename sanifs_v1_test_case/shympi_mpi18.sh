#!/bin/bash

#BSUB -x 
#BSUB -q p_short
#BSUB -J sani64
#BSUB -n 18
#BSUB -o LOGS/logout.%J.out
#BSUB -e LOGS/logerr.%J.err
#BSUB -P R000

module purge
module purge
module load intel19.5/19.5.281
module load impi19.5/19.5.281
module load impi19.5/petsc/3.7.5
export RUNTIME_OPTS="-ksp_type bcgs -ksp_rtol 1e-18 -ksp_atol 1e-15  "

if [ "$I_MPI_HYDRA_BOOTSTRAP" == "" ]; then
  export I_MPI_HYDRA_BOOTSTRAP=lsf
fi
export I_MPI_HYDRA_BRANCH_COUNT=1

if [ "$I_MPI_HYDRA_BRANCH_COUNT" == "" ]; then
  export I_MPI_HYDRA_BRANCH_COUNT=`cat $LSB_DJOB_HOSTFILE | uniq | wc -l`
fi
if [ "$I_MPI_LSF_USE_COLLECTIVE_LAUNCH" == "" ]; then
  export I_MPI_LSF_USE_COLLECTIVE_LAUNCH=1
fi

mpiexec.hydra -l /users_home/opa/ses-dev/shyfem_versions/shympi_r8_4/fem3d/shympi param_mpi18.str $RUNTIME_OPTS

echo "Job completed at: " `date`

sleep 10 

