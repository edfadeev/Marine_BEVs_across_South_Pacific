#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
#SBATCH --mail-user=ALL
#SBATCH --output=./00_LOGS/%x-%j.out

#export functions of each gene
anvi-export-functions -c $WORKDIR/03_CONTIGS/$PROJECT-contigs.db --annotation-sources $db \
-o $WORKDIR/03_CONTIGS/$PROJECT-$db-functions.txt

#remove spaces
sed -i -e 's/ /_/g' $WORKDIR/03_CONTIGS/$PROJECT-$db-functions.txt 

#subset to only relevant gene calls
awk -F' ' 'NR==FNR{gcids[$1];next} $1 in gcids' $WORKDIR/11_PROTEIN/R_data/prot_GCIDs.txt $WORKDIR/03_CONTIGS/$PROJECT-$db-functions.txt > $WORKDIR/11_PROTEIN/$PROJECT-detected-proteins_$db.txt

rm $WORKDIR/03_CONTIGS/$PROJECT-$db-functions.txt