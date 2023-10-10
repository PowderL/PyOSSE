#!/bin/bash

MEM=`ulimit -s`

if [ $MEM == 10240 ]; then
    MEM=100000000
fi
MEM=999999999
# setting to unlimited means that the code doesn't work.
ulimit -s $MEM
MEM=`ulimit -s`

# echo $MEM
# export STACKSIZE=$MEM
export KMP_STACKSIZE=$MEM
export OMP_NUM_THREADS=4
export NCPUS=4

# list=`sar -P ALL 10 | grep Average |cut -c26-28` 
# list=`ls -l *.h |cut -c20-22`

usage1=0
usage2=0
usage3=0
usage4=0
i=-1

echo 'run with flag ', $run_flag

# run_flag="--interleave 3  --cpunodebind 3"

# numactl $run_flag  ./geos_3h_nep

# numactl $run_flag  ./geos_3h

./geos_ncep


# numactl $run_flag  ./geos_daily

# numactl $run_flag  ./geos_dd_bb

# numactl $run_flag  ./geos_mm_bb











