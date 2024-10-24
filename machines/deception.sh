# This file loads all necessary modules for building Haero on Deception (PNNL)
. /etc/profile.d/modules.sh
module purge

module load gcc/11.2.0
module load cmake/3.28.1
module load python/3.7.0
module load cuda/11.7
