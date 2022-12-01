# This file loads all necessary modules for building Haero on Deception (PNNL)
. /etc/profile.d/modules.sh
module purge 

module load gcc/9.1.0
module load cmake/3.21.4
module load python/3.7.0
module load cuda/11.4
