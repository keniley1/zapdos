#!/bin/bash
#PBS -M keniley1@illinois.edu
#PBS -m abe
#PBS -N water_ph3_003
#PBS -P moose
#PBS -l select=1:ncpus=12:mpiprocs=12
#PBS -l walltime=8:00:00

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE

#module load use.moose PETSc
module load use.moose PETSc/3.11.4-GCC
MV2_ENABLE_AFFINITY=0 mpiexec ~/projects/lemhi/zapdos/zapdos-opt -i sim_ph_003.i --color off >> sim_ph3_003_stdout
