#!/bin/bash
#SBATCH -p batch
#SBATCH -N 20
#SBATCH -n 20
#SBATCH --time=08:30:00
#SBATCH --mem=1G

# Notification configuration
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a1816653@adelaide.edu.au

#module load Python/3.9.6-GCCcore-11.2.0
#module load SciPy-bundle/2021.10-foss-2021b
#module load matplotlib/3.5.2-foss-2021b

module load Anaconda3/2023.03

# Activate the conda environment
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate myenv
python -c "import shapely; print('Shapely version:', shapely.__version__)"

python run_pso_mcts.py
