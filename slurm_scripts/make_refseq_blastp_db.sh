#!/bin/bash
#
#SBATCH --cpus-per-task=2
#SBATCH --mem=100GB
#SBATCH --partition=basic
#SBATCH --mail-user=eduard.fadeev@univie.ac.at
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=4:00:00

#load module
module load ncbiblastplus/2.16.0

# Directory that will contain the blastp index
DBDIR=/lisc/scratch/oceanography/efadeev/resources/ref_dbs/NCBI

cd $TMPDIR

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary_refseq.txt > ftpdirpaths
awk 'BEGIN{FS=OFS="/";filesuffix="protein.faa.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
wget -i ftpfilepaths
gunzip *.gz
cat *.faa > bac_refseq.fa
makeblastdb -in bac_refseq.fa -out bacterial_refseq -dbtype prot


mkdir $DBDIR/bac_refseq
mv $TMPDIR/bacterial_refseq* $DBDIR/bac_refseq/

rm -Rf $TMPDIR