#!/bin/bash
#
#SBATCH --cpus-per-task=2
#SBATCH --mem=10GB
#SBATCH --partition=basic
#SBATCH --mail-user=eduard.fadeev@univie.ac.at
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=1:00:00

#load module
module load BLAST+/2.17.0-gompi-2025a

#generate array with protein chunk names
i=0
while read line
do
chunks[$i]="$line"
i=$((i+1))
done < $OUTDIR/chunk.list


# Exit the slurm script if a command fails
set -e

# Name of the blastp index
DBNAME=TCDB_ref

for chunk in ${chunks[${SLURM_ARRAY_TASK_ID}]};
do
file=$(basename $chunk)

cp $chunk $TMPDIR/input

psiblast -db $OUTDIR/$DBNAME/$DBNAME -query $TMPDIR/input -max_target_seqs 10 \
-outfmt "6 qseqid sseqid pident length evalue bitscore" -out $TMPDIR/blast.out -num_threads $SLURM_CPUS_PER_TASK

cp $TMPDIR/blast.out $OUTDIR_TCDB/${file}_${DBNAME}.out

done

rm -rf $TMPDIR
