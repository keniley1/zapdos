#!/bin/bash
#PBS -M keniley1@illinois.edu
#PBS -m abe
#PBS -N Ar_050_humid_02
#PBS -P moose
#PBS -l select=1:ncpus=12:mpiprocs=12
#PBS -l walltime=08:00:00

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE

#module load use.moose PETSc
module load use.moose PETSc/3.11.4-GCC
MV2_ENABLE_AFFINITY=0 mpiexec ~/projects/lemhi/zapdos/zapdos-opt -i sim_h050.i --color off >> stdout_h050_02
