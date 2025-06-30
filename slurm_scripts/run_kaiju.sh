#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --mem=120GB
#SBATCH --partition=basic
#SBATCH --mail-user=ALL
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=4:00:00

# Exit the slurm script if a command fails
set -e

# Directory that contains the blastp index. Replace "2025_01" with recent release
DBDIR=/lisc/scratch/oceanography/efadeev/anvio-resources/kaiju_dbs/

# Name of the blastp index
DBNAME=kaiju_db_refseq_2024-08-14

# We call lisc_localcache with DBDIR as argument
# The function returns the path of the local cache - we store it in the DBCACHE variable 
DBCACHE=$(lisc_localcache $DBDIR)

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

#add gene taxonomic assignments to genes
kaiju -t $DBDIR/$DBNAME/nodes.dmp \
-f $DBDIR/$DBNAME/kaiju_db_refseq.fmi \
-i $chunk \
-o $TMPDIR/gene_calls_refseq_tax.out \
-z $SLURM_CPUS_PER_TASK \
-v

#add taxon name
kaiju-addTaxonNames -t $DBDIR/$DBNAME/nodes.dmp \
-n $DBDIR/$DBNAME/names.dmp \
-i $TMPDIR/gene_calls_refseq_tax.out \
-o $TMPDIR/gene_calls_refseq_names.out \
-r superkingdom,phylum,order,class,family,genus,species

mv $TMPDIR/gene_calls_refseq_tax.out $WORKDIR/11_PROTEIN/temp_taxa/${file}_refseq_tax.out

mv $TMPDIR/gene_calls_refseq_names.out $WORKDIR/11_PROTEIN/temp_taxa/${file}_refseq_names.out

done 

rm -Rf $TMPDIR
