#!/bin/bash
#
#SBATCH --cpus-per-task=16
#SBATCH --mem=120GB
#SBATCH --partition=basic
#SBATCH --mail-user=eduard.fadeev@univie.ac.at
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=12:00:00


CAT_PATH=/lisc/data/scratch/oceanography/efadeev/resources/CAT_pack/CAT_pack
DB_PATH=/lisc/data/scratch/oceanography/efadeev/resources/ref_dbs/CAT_dbs/20241212_CAT_nr_website/db
TAX_PATH=/lisc/data/scratch/oceanography/efadeev/resources/ref_dbs/CAT_dbs/20241212_CAT_nr_website/tax

cp $WORKDIR/13_PD_masked/contigs_of_interest.fa $TMPDIR/contigs
cp $WORKDIR/13_PD_masked/proteins_in_contigs.faa $TMPDIR/proteins

pushd $TMPDIR

$CAT_PATH/CAT_pack contigs --nproc $SLURM_CPUS_PER_TASK --block_size 8 \
-c $TMPDIR/contigs \
-d $DB_PATH -t $TAX_PATH -p $TMPDIR/proteins -o $TMPDIR/CAT_output

cp $TMPDIR/CAT_output* $WORKDIR/13_PD_masked/

popd
