#!/bin/bash
#
#SBATCH --cpus-per-task=2
#SBATCH --mem=10GB
#SBATCH --partition=basic
#SBATCH --mail-user=eduard.fadeev@univie.ac.at
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=2:00:00

#load module
module load java/11.0.4

export PATH=/lisc/data/scratch/oceanography/efadeev/resources/interproscan-5.66-98.0:$PATH

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

interproscan.sh -i $chunk -f tsv --cpu $SLURM_CPUS_PER_TASK \
--applications $db --disable-precalc --iprlookup --goterms -o $TMPDIR/interpro-output.tsv

mv $TMPDIR/interpro-output.tsv $OUTDIR_INTER/${file}_${db}_interpro-output.tsv

done
