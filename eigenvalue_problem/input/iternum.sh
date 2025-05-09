#!/bin/zsh

mpirun -np 16 ./iter_num \
  -ksp_type cg \
  -pc_type none \