#!/bin/bash
#
#SBATCH --cpus-per-task=32
#SBATCH --mem=120GB
#SBATCH --partition=basic
#SBATCH --mail-user=ALL
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=24:00:00

# Exit the slurm script if a command fails
set -e

# Directory that contains the blastp index. Replace "2025_01" with recent release
DBDIR=/lisc/scratch/oceanography/efadeev/anvio-resources/kaiju_dbs/

# Name of the blastp index
DBNAME=kaiju_db_refseq_2024-08-14

# We call lisc_localcache with DBDIR as argument
# The function returns the path of the local cache - we store it in the DBCACHE variable 
DBCACHE=$(lisc_localcache $DBDIR)


R1=$(ls -1 ./01_QC/*R1.fastq.gz| awk '{printf "%s,", $0}' - | sed 's/,$//')
R2=$(ls -1 ./01_QC/*R2.fastq.gz| awk '{printf "%s,", $0}' - | sed 's/,$//')

#add gene taxonomic assignments to genes
kaiju-multi -t $DBDIR/$DBNAME/nodes.dmp \
-f $DBDIR/$DBNAME/kaiju_db_refseq.fmi \
-i $R1 \
-j $R2 \
-z $SLURM_CPUS_PER_TASK > $TMPDIR/fastq_refseq_tax.out \

#add taxon name
kaiju-addTaxonNames -t $DBDIR/$DBNAME/nodes.dmp \
-n $DBDIR/$DBNAME/names.dmp \
-i $TMPDIR/fastq_refseq_tax.out \
-o $TMPDIR/fastq_refseq_names.out \
-r superkingdom,phylum,order,class,family,genus,species

mv $TMPDIR/fastq_refseq_tax.out $WORKDIR/05_TAXONOMY/fastq_refseq_tax.out

mv $TMPDIR/fastq_refseq_names.out $WORKDIR/05_TAXONOMY/fastq_refseq_names.out

rm -Rf $TMPDIR
