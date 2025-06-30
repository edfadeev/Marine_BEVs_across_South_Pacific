#!/bin/bash
#
#SBATCH --cpus-per-task=4
#SBATCH --mem=20GB
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
DBDIR=/lisc/scratch/oceanography/efadeev/ref_dbs/merops

# Name of the blastp index
DBNAME=pepunit

# We call lisc_localcache with DBDIR as argument
# The function returns the path of the local cache - we store it in the DBCACHE variable 
DBCACHE=$(lisc_localcache $DBDIR)

for chunk in ${chunks[${SLURM_ARRAY_TASK_ID}]};
do
file=$(basename $chunk)

cp $chunk $TMPDIR/input

blastp -db $DBCACHE/$DBNAME -query $TMPDIR/input -max_target_seqs 1 \
-outfmt "6 qseqid sseqid pident length evalue bitscore" -out $TMPDIR/blastp.out -num_threads $SLURM_CPUS_PER_TASK

#filter based on identity of 30%, evalue < 0.001 and bit score >50 and then merge with annotations
awk '$3 > 25 && $5 < 0.001 && $6 > 50 {print $1,$2}' $TMPDIR/blastp.out | \
awk 'NR == FNR {x[$1]=$2; next}{ print $0,x[$2]}' \
/lisc/scratch/oceanography/efadeev/ref_dbs/merops/pepunit.headers - |\
sed 's/\s*$/ MEROPS/' - |  awk 'BEGIN {FS=" "; OFS="\t"} {print $1,$5,$2,$4,$3}' - > $TMPDIR/blastp.ann

mv $TMPDIR/blastp.ann $WORKDIR/11_PROTEIN/temp_merops/${file}_blastp_merops.ann

done

rm -rf $TMPDIR