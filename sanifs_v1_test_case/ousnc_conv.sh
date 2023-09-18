#!/bin/bash
#BSUB -q s_short 
#BSUB -J dm_ousnc
#BSUB -o conv_logs/ousnc.%J_%I.out
#BSUB -e conv_logs/ousnc.%J_%I.err
#BSUB -P R000

module purge
module purge
module load intel19.5/19.5.281
module load impi19.5/19.5.281
module load impi19.5/petsc/3.7.5
module load intel19.5/netcdf/C_4.7.2-F_4.5.2_CXX_4.3.1

sim_ous_name="sani64_mpi1_semi_implicit_final_64bit.ous"
sim_bas_name="newsani64.bas"

output_name=`echo $sim_ous_name | rev | cut -c 5- | rev`

echo "output name=" $output_name

memory -s $sim_ous_name
memory -b $sim_bas_name

apn2nc="apnous2nc.str"

rm $apn2nc

echo -e " " >> "$apn2nc"
echo -e " " >> "$apn2nc"
echo -e " " >> "$apn2nc" #window tot
echo -e " " >> "$apn2nc"
echo -e " " >> "$apn2nc"
echo -e " " >> "$apn2nc"

mpiexec.hydra -l /users_home/opa/ses-dev/shyfem_versions/2nc/ous2nc < "$apn2nc"

sleep 5

mv out.nc ${output_name}.ous.nc
