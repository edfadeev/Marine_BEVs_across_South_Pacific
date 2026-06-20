#!/bin/bash
#
#SBATCH --cpus-per-task=6
#SBATCH --mem=3GB
#SBATCH --mail-user=eduard.fadeev@univie.ac.at
#SBATCH --partition=basic
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=00:30:00

#load module 
module load SignalP/6.0i-foss-2024a-fast

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

signalp6 --fastafile $chunk --output_dir $TMPDIR/${file} --format txt --torch_num_threads $SLURM_CPUS_PER_TASK

mv $TMPDIR/${file} $OUTDIR_SIGNALP/

done
