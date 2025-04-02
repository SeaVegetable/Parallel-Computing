#!/bin/zsh

mpirun -np 9 ./iter_num_fem \
  -ksp_type gmres \
  -pc_type hypre \
  -ksp_rtol 1.0e-6 \
  -ksp_atol 1.0e-50 \
  -ksp_max_it 200 \
  -ksp_gmres_restart 200 \
  -pc_hypre_boomeramg_coarsen_type HMIS \
  -pc_hypre_boomeramg_interp_type ext+i \
  -pc_hypre_boomeramg_truncfactor 0.3 \
  -pc_hypre_boomeramg_strong_threshold 0.5 \
  -pc_hypre_boomeramg_P_max 5 \
  -pc_hypre_boomeramg_agg_nl 2 \
  -log_view > iternum_fem.log 2>&1 \