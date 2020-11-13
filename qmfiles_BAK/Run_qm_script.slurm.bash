#!/bin/bash
##SBATCH -D  /home/yichun/cx26_K_milestoning/test_amber20/
#SBATCH -J psi4QM
#SBATCH --partition=defq
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
##SBATCH --gres=gpu:1
#SBATCH --time=08:00:00
#####SBATCH --share

eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda init bash
source activate psi4
which python

scriptName=$1
logFileName=`echo $1 | sed "s:.py:.out:g"`

echo "Running $scriptName. Saving result to $logFileName"
python $scriptName &> $logFileName

echo "Done"
