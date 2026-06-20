#!/bin/bash
#
#SBATCH --cpus-per-task=2
#SBATCH --mem=10GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --partition=basic
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=8:00:00

#load module 
#module unload python3
#module load deeplocpro

module load Conda; conda activate DeepLocPro-1.0.0

# Exit the slurm script if a command fails
set -e

#generate array with protein chunk names
i=0
while read line
do
chunks[$i]="$line"
i=$((i+1))
done < $OUTDIR/chunk.list


for chunk in ${chunks[${SLURM_ARRAY_TASK_ID}]};
do
file=$(basename $chunk)

deeplocpro -f $chunk -o $TMPDIR/${file}

mv $TMPDIR/${file} $OUTDIR_DEEPLOC/

done
