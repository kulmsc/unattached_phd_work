#! /bin/bash -l

#SBATCH --partition=panda_physbio   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=rgwas
#SBATCH --time=72:00:00   # HH/MM/SS
#SBATCH --mem=24G   # memory requested, units available: K,M,G,T

source ~/.bashrc

spack load -r /cypfv4n miniconda3@4.3.14%gcc@6.3.0 arch=linux-centos7-x86_64

source activate basic

python general_prev.py 350000 400000


