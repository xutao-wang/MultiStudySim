#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xutaowang@hsph.harvard.edu
#SBATCH -J base_5s_10p
#SBATCH -n 18  
#SBATCH -N 1    
#SBATCH -p general  
#SBATCH --mem 32G
#SBATCH -t 0-12:00
#SBATCH -o test-%x.%j-%a.out 
#SBATCH -e test-%x.%j-%a.err 

module load R/3.5.1-fasrc01
export R_LIBS_USER=$HOME/apps/R_3.5.1:$R_LIBS_USER

Rscript --vanilla ~/MultiSim/result/beta_sigma_5s_10p/Interaction/oneInt/beta_sigma_5s_10p_oneInt.R
