#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

in=0
let "out=in+1"

time ../runmd.py \
     --xml system.xml --restrain ':1-129' \
     --temp 100.0  --nrespa 2 --dt 2.0 \
     -p 4lzt_uc.pdb -x md$out.nc -r md$out.x -o md$out.o \
     -s md$in.x -n 10000 --interval 500 --gamma_ln 5.0 > md$out.log 2>&1 

