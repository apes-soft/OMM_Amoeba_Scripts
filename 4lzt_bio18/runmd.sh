#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

in=0
let "out=in+1"

../runmd.py \
     --xml system.xml --restrain ':1-129,414-419' -k 10.0 \
     --temp 100.0  --nrespa 2 --dt 1.0 \
     -p 4lzt_uc.pdb_2 -x md$out.nc -r md$out.x -o md$out.o \
     -s md$in.x -n 100000 --interval 10000 --gamma_ln 5.0 > md$out.log 2>&1 

