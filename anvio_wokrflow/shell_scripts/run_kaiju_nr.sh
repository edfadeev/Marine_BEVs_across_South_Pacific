#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --mem=220GB
#SBATCH --partition=basic
#SBATCH --mail-user=ALL
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=4:00:00

# Exit the slurm script if a command fails
set -e

# Directory that contains the blastp index. Replace "2025_01" with recent release
DBDIR=/lisc/data/scratch/oceanography/efadeev/resources/ref_dbs/kaiju_dbs/

# Name of the blastp index
#DBNAME=kaiju_db_refseq_2024-08-14
DBNAME=kaiju_db_nr_2024-08-25

# We call lisc_localcache with DBDIR as argument
# The function returns the path of the local cache - we store it in the DBCACHE variable 
DBCACHE=$(lisc_localcache $DBDIR)

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

#add gene taxonomic assignments to genes
kaiju -t $DBDIR/$DBNAME/nodes.dmp \
-f $DBDIR/$DBNAME/kaiju_db_nr.fmi \
-i $chunk \
-o $TMPDIR/gene_calls_${DBNAME}_tax.out \
-z $SLURM_CPUS_PER_TASK \
-v -p

#add taxon name
kaiju-addTaxonNames -t $DBDIR/$DBNAME/nodes.dmp \
-n $DBDIR/$DBNAME/names.dmp \
-i $TMPDIR/gene_calls_${DBNAME}_tax.out \
-o $TMPDIR/gene_calls_${DBNAME}_names.out \
-r superkingdom,phylum,order,class,family,genus,species

mv $TMPDIR/gene_calls_${DBNAME}_tax.out $OUTDIR_TAXA/${file}_${DBNAME}_tax.out

mv $TMPDIR/gene_calls_${DBNAME}_names.out $OUTDIR_TAXA/${file}_${DBNAME}_names.out

done 
