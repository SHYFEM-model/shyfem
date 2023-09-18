#!/bin/bash
#BSUB -q s_short
#BSUB -J dm_nosnc
#BSUB -o conv_logs/nosnc.%J_%I.out
#BSUB -e conv_logs/nosnc.%J_%I.err
#BSUB -P R000

module purge
module purge
module load intel19.5/19.5.281
module load impi19.5/19.5.281
module load impi19.5/petsc/3.7.5
module load intel19.5/netcdf/C_4.7.2-F_4.5.2_CXX_4.3.1

sim_nos_name="sani64_mpi1_semi_implicit_final_64bit.nos"
sim_bas_name="newsani64.bas"

output_name=`echo $sim_nos_name | rev | cut -c 5- | rev`

echo "output name="$output_name

memory -s $sim_nos_name
memory -b $sim_bas_name

apn2nc="apnnos2nc.str"

rm $apn2nc

echo -e " " >> "$apn2nc"
echo -e " " >> "$apn2nc"
echo -e " " >> "$apn2nc"
echo -e " " >> "$apn2nc"
echo -e " " >> "$apn2nc"
echo -e " " >> "$apn2nc"

mpiexec.hydra -l /users_home/opa/ses-dev/shyfem_versions/2nc/nos2nc < "$apn2nc"

sleep 5

mv out.nc ${output_name}.nos.nc
