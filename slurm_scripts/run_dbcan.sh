#!/bin/bash
#
#SBATCH --cpus-per-task=4
#SBATCH --mem=20GB
#SBATCH --partition=basic
#SBATCH --mail-user=eduard.fadeev@univie.ac.at
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=4:00:00

#load module
module load dbcan/5.1.1-3.13.3

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
DBDIR=/lisc/scratch/mirror/dbcan

# Name of the blastp index
DBNAME=5.0.6

# We call lisc_localcache with DBDIR as argument
# The function returns the path of the local cache - we store it in the DBCACHE variable 
DBCACHE=$(lisc_localcache $DBDIR)

for chunk in ${chunks[${SLURM_ARRAY_TASK_ID}]};
do
file=$(basename $chunk)

run_dbcan CAZyme_annotation --input_raw_data $chunk --mode protein --output_dir $TMPDIR/${file} \
--db_dir $DBCACHE/$DBNAME --threads $SLURM_CPUS_PER_TASK

awk 'NR > 1 && $6 > 1 {print $1, $2, $7}' $TMPDIR/${file}/overview.tsv > \
$TMPDIR/${file}_dbcan_results.tsv

mv $TMPDIR/${file}_dbcan_results.tsv $WORKDIR/11_PROTEIN/temp_dbcan/${file}_dbcan_results.tsv

done

rm -rf $TMPDIR