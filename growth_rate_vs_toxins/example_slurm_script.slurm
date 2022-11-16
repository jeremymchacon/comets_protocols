#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=10g
#SBATCH --mail-type=NONE

module load conda
module load java/openjdk-13.0.2
module load gurobi/9.0.2
unset PYTHONPATH
unset PYTHONHOME
unset PYTHONSTARTUP

source activate streptomyces

cd RPS

python using_python_scripts_from_command_line_with_arguments.py $SLURM_ARRAY_TASK_ID

