#!/bin/bash
#  run_corr_ind.sh
#  
#
#  Created by A Lyne on 14/05/2018.
#
#PBS -l nodes=1:ppn=1

# set max wallclock time
#PBS -l walltime=0:30:00

# set memory requirements
#PBS -l mem=10gb

# set name of job
#PBS -N tree_comp_res

#PBS -k oe

R_EXEC="/bioinfo/local/build/Centos/R/R-3.5.0/bin/Rscript"
R_FILE="${PBS_O_WORKDIR}/compile_results.R"
$R_EXEC $R_FILE
