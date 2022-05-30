#!/bin/bash
#PBS -q gold5120
#PBS -l nodes=1:ppn=1
#PBS -o job.log
#PBS -e job.err
ulimit -s unlimited
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1
EXEC=minimax_grid_main.py
python3 $EXEC
