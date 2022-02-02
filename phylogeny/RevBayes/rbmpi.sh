#!/bin/bash
#SBATCH --account=adamians   ### change this to your actual account for charging
#SBATCH --partition=short       ### queue to submit to
#SBATCH --job-name=GVrbTopology.job    ### job name
#SBATCH --output=GVrbTopology.out   ### file in which to store job stdout
#SBATCH --error=GVrbTopology.err    ### file in which to store job stderr
#SBATCH --time=3-00:00:00              ### wall-clock time limit
#SBATCH --mem=1G             ### memory limit per node, in MB
#SBATCH --nodes=2               ### number of nodes to use
#SBATCH --ntasks-per-node=2     ### number of tasks to launch per node
#SBATCH --cpus-per-task=2       ### number of cores for each task