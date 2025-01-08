#!/bin/bash
#SBATCH --time=21:00:00
#SBATCH --partition=single
#SBATCH --nodes=1
#SBATCH --cpus-per-task=80
#SBATCH --mem=20000mb
#SBATCH --export=ALL
#SBATCH --output=sim_results/MCCW.out
#SBATCH --error=sim_results/MCCW.err
#SBATCH -J gprMax_MCCW_plate_no_obstacle

# Add the following commands for email notifications:
# #SBATCH --mail-type=all
# #SBATCH --mail-user=your@email-adress

# Directory of gprMax, adjust according to your folder structure
cd ..
cd ..
cd ..

module load compiler/gnu/10.2
module load devel/miniconda/4.9.2

#Usually you should set
export KMP_AFFINITY=compact,1,0
#export KMP_AFFINITY=verbose,compact,1,0 prints messages concerning the supported affinity
#KMP_AFFINITY Description: https://software.intel.com/en-us/node/524790#KMP_AFFINITY_ENVIRONMENT_VARIABLE

export OMP_NUM_THREADS=$((${SLURM_JOB_CPUS_PER_NODE}/2))

# Adjust parameter after -n according to the amount of performed simulations
startsim="python -m gprMax bioradar_models/Monochromatic_Continuous_Wave/Sim2_plate_no_obstacle/MCCW_radar.in -n 101"

conda init bash
source ~/.bashrc
conda activate gprMax
exec $startsim
conda deactivate