#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --mem=200GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --partition=basic
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=12:00:00

#cluster proteins using CD-HIT
module load cdhit/4.8.1

anvi-get-sequences-for-gene-calls --contigs-db $WORKDIR/03_CONTIGS/$PROJECT-contigs.db \
--get-aa-sequences --output-file $TMPDIR/proteins.fasta

cd-hit -i $TMPDIR/proteins.fasta \
-aS 0.9 -c 0.99 -o $TMPDIR/$PROJECT-proteins_clust_99.fasta \
-M $SLURM_MEM_PER_NODE -T $SLURM_CPUS_PER_TASK -G 0

cd-hit -i $TMPDIR/proteins.fasta \
-aS 0.9 -c 0.95 -o $TMPDIR/$PROJECT-proteins_clust_95.fasta \
-M $SLURM_MEM_PER_NODE -T $SLURM_CPUS_PER_TASK -G 0

mv $TMPDIR/$PROJECT-proteins_clust_99.fasta $TMPDIR/$PROJECT-proteins_clust_95.fasta \
$WORKDIR/11_PROTEIN/

rm -Rf $TMPDIR