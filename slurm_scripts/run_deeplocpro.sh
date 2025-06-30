#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --partition=basic
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=2:00:00

#load module 
module load deeplocpro

# Exit the slurm script if a command fails
set -e

#generate array with protein chunk names
i=0
while read line
do
chunks[$i]="$line"
i=$((i+1))
done < $WORKDIR/11_PROTEIN/chunk.list


for chunk in ${chunks[${SLURM_ARRAY_TASK_ID}]};
do
file=$(basename $chunk)

deeplocpro -f $chunk -d cuda -o $TMPDIR/${file}

mv $TMPDIR/${file} $WORKDIR/11_PROTEIN/temp_deeploc/

done

rm -rf $TMPDIR