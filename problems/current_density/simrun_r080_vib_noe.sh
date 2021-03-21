#!/bin/bash
#PBS -M keniley1@illinois.edu
#PBS -m abe
#PBS -N sim_r080_vib_noe_01
#PBS -P moose
#PBS -l select=1:ncpus=12:mpiprocs=12
#PBS -l walltime=08:00:00

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE

#module load use.moose PETSc
module load use.moose PETSc/3.11.4-GCC
MV2_ENABLE_AFFINITY=0 mpiexec ~/projects/lemhi/zapdos/zapdos-opt -i sim_r080_vib_noe.i --color off >> stdout_r080_vib_noe_01
