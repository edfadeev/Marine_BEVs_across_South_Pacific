#!/bin/bash
#
#SBATCH --cpus-per-task=16
#SBATCH --mem=70GB
#SBATCH --partition=basic
#SBATCH --mail-user=eduard.fadeev@univie.ac.at
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=18:00:00

#load module
#module load ncbiblastplus/2.16.0
module load BLAST+/2.17.0-gompi-2025a


#add annotation parser to path
export PATH=/lisc/data/scratch/oceanography/efadeev/resources/bioinf_scripts:$PATH

#generate array with protein chunk names
i=0
while read line
do
chunks[$i]="$line"
i=$((i+1))
done < $OUTDIR/chunk.list


# Exit the slurm script if a command fails
set -e

# Directory that contains the blastp index. Replace "2025_01" with recent release
DBDIR=/lisc/data/scratch/oceanography/efadeev/resources/ref_dbs/NCBI/2025_05_12/

# Name of the blastp index
DBNAME=bacterial_refseq

# We call lisc_localcache with DBDIR as argument
# The function returns the path of the local cache - we store it in the DBCACHE variable 
DBCACHE=$(lisc_localcache $DBDIR)

for chunk in ${chunks[${SLURM_ARRAY_TASK_ID}]};
do
file=$(basename $chunk)

cp $chunk $TMPDIR/input

blastp -db $DBCACHE/$DBNAME -query $TMPDIR/input -max_target_seqs 10 \
-outfmt "6 qseqid sseqid pident length evalue bitscore" -out $TMPDIR/blastp.out -num_threads $SLURM_CPUS_PER_TASK

cp $TMPDIR/blastp.out $OUTDIR_BLASTP/${file}_${DBNAME}_blastp_mult.out

annotate_blast_hits $TMPDIR/blastp.out p eduard.fadeev@univie.ac.at

mv $TMPDIR/blastp.out_annotated.tsv $OUTDIR_BLASTP/${file}_${DBNAME}_annotated_mult.tsv

done
