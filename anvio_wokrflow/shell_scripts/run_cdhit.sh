#!/bin/bash
#
#SBATCH --cpus-per-task=16
#SBATCH --mem=200GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --partition=basic
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=12:00:00

#cluster proteins using CD-HIT
module load CD-HIT

anvi-get-sequences-for-gene-calls --contigs-db $WORKDIR/03_CONTIGS/$PROJECT-contigs.db \
--get-aa-sequences --output-file $TMPDIR/$PROJECT-proteins.fasta

python $REPO_DIR/anvio_workflow/shell_scripts/filter_tryptic.py $TMPDIR/$PROJECT-proteins.fasta $TMPDIR/$PROJECT-proteins_tryptic.fasta

cd-hit -i $TMPDIR/$PROJECT-proteins_tryptic.fasta \
-g 1 -c 1.0 -o $TMPDIR/$PROJECT-proteins_nr100.fasta \
-M $SLURM_MEM_PER_NODE -T $SLURM_CPUS_PER_TASK

mv $TMPDIR/* $OUTDIR/

rm -Rf $TMPDIR
