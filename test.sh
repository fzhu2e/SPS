#!/bin/bash
./makemake.perl sp_dynamic
make clean; make
ulimit -s unlimited
#OMP_NUM_THREADS=1
#OMP_NUM_THREADS=2
OMP_NUM_THREADS=4
#OMP_NUM_THREADS=8
time ./sp_dynamic
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"

exit
