#! /bin/bash -l
 
#SBATCH --partition=panda_physbio   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=rgwas
#SBATCH --time=72:00:00   # HH/MM/SS
#SBATCH --mem=14G   # memory requested, units available: K,M,G,T
 
source ~/.bashrc
spack load -r /mjrrusu 

Rscript residual_compare.R $1


