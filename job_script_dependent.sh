#!/bin/sh

#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=240G
#SBATCH -o log/%x_%j.out
#SBATCH --mail-user=frederic.ouimet.23@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=message/%x_%j-result.txt
#SBATCH --error=message/%x_%j-error.txt
#SBATCH --exclusive
#SBATCH --account=def-fouimet

module load StdEnv/2023 r/4.3.1
Rscript simulations_dep_parallel.R
