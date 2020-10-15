#!/bin/bash
#PBS -M keniley1@illinois.edu
#PBS -m abe
#PBS -N air_R01MOhm_01
#PBS -P moose
#PBS -l select=1:ncpus=12:mpiprocs=12:mem=40gb
#PBS -l walltime=48:00:00

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE

#module load use.moose PETSc
module load use.moose PETSc/3.11.4-GCC
MV2_ENABLE_AFFINITY=0 mpiexec ~/projects/zapdos/zapdos-opt -i air_water_R01MOhm.i --color off >> air_water_R01MOhm_01_out
