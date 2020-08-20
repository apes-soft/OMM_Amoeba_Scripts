#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

in=0
let "out=in+1"

time  ../runmd.py \
           --xml system.xml \
           -p 2igd_solvated.pdb -x 2igd.md$out.nc \
           -r 2igd.md$out.x -o 2igd.md$out.o \
           -s 2igd.md$in.x -n 10000 --interval 500 \
           --gamma_ln 5.0 > 2igd.md$out.log 2>&1 

