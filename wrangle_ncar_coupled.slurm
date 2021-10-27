#!/bin/bash
#
#SBATCH --job-name=wrangle_coupled_data
#SBATCH --output=tempout_coupled.txt
#SBATCH --ntasks-per-node=28
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH -p long-28core

module load shared
module load netcdf/gcc/4.7.1
module load gdal/2.2.3
module load R/3.6.2
module load gcc-stack
module load openblas/dynamic/0.2.18
module load lapack/gcc/64/3.6.0
export R_LIBS=/gpfs/home/bsen/R_libs

cd /gpfs/home/bsen/NCAR_data_wrangling
Rscript wrangle_ncar_coupled.R