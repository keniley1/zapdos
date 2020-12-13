#!/bin/bash
#PBS -M keniley1@illinois.edu
#PBS -m abe
#PBS -N water_2d_001
#PBS -P moose
#PBS -l select=1:ncpus=24:mpiprocs=24
#PBS -l walltime=48:00:00

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE

#module load use.moose PETSc
module load use.moose PETSc/3.11.4-GCC
MV2_ENABLE_AFFINITY=0 mpiexec ~/projects/lemhi/zapdos/zapdos-opt -i sim2d.i --color off >> water_2d_stdout_001
