#!/bin/bash
#PBS -M keniley1@illinois.edu
#PBS -m abe
#PBS -N h2o2_neutral_0mM_01
#PBS -P moose
#PBS -l select=1:ncpus=12:mpiprocs=12
#PBS -l walltime=24:00:00

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE

#module load use.moose PETSc
module load use.moose PETSc/3.11.4-GCC
MV2_ENABLE_AFFINITY=0 mpiexec ~/projects/lemhi/zapdos/zapdos-opt -i sim_h2o2_neutral_0mM.i --color off >> h2o2_neutral_0mM_terminal_01
