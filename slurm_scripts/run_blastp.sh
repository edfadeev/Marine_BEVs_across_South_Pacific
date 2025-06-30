#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --mem=320GB
#SBATCH --partition=basic
#SBATCH --mail-user=eduard.fadeev@univie.ac.at
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=4:00:00

#load module
module load ncbiblastplus/2.16.0

#add annotation parser to path
export PATH=/lisc/scratch/oceanography/efadeev/anvio-resources/bioinf_scripts:$PATH

#generate array with protein chunk names
i=0
while read line
do
chunks[$i]="$line"
i=$((i+1))
done < $WORKDIR/11_PROTEIN/chunk.list


# Exit the slurm script if a command fails
set -e

# Directory that contains the blastp index. Replace "2025_01" with recent release
DBDIR=/lisc/scratch/mirror/ncbi/2025-03-27

# Name of the blastp index
DBNAME=nr

# We call lisc_localcache with DBDIR as argument
# The function returns the path of the local cache - we store it in the DBCACHE variable 
DBCACHE=$(lisc_localcache $DBDIR)

for chunk in ${chunks[${SLURM_ARRAY_TASK_ID}]};
do
file=$(basename $chunk)

cp $chunk $TMPDIR/input

blastp -db $DBCACHE/$DBNAME -query $TMPDIR/input -max_target_seqs 1 \
-outfmt "6 qseqid sseqid pident length evalue bitscore" -out $TMPDIR/blastp.out -num_threads $SLURM_CPUS_PER_TASK

cp $TMPDIR/blastp.out $WORKDIR/11_PROTEIN/temp_blastp/${file}_blastp.out

annotate_blast_hits $TMPDIR/blastp.out p eduard.fadeev@univie.ac.at

mv $TMPDIR/blastp.out_annotated.csv $WORKDIR/11_PROTEIN/temp_blastp/${file}_blastp_annotated.csv

done

rm -rf $TMPDIR