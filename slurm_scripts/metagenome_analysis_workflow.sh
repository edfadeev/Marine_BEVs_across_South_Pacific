# Anvio workflow for metagenomes assembly, annotation and binning
#login with port forwarding: ssh -L 5678:localhost:5678 slurm

################################
# Define variables
################################
export PATH=/home/user/fadeev/anvio-resources/centrifuge_tool/bin:/home/user/fadeev/anvio-resources/DAS_Tool/:/home/user/fadeev/anvio-resources/pullseq/src/:/lisc/user/fadeev/anvio-resources/InterProScanParser:/home/user/fadeev/anvio-resources/KrakenTools:/lisc/user/fadeev/.local/bin:$PATH

#Define working directory
WORKDIR=/lisc/scratch/oceanography/efadeev/SO289/SO289_coassembly/ && cd $WORKDIR
PROJECT=SO289
REPO_DIR=/lisc/scratch/oceanography/efadeev/anvio_scripts

################################################################
# Run metagenome assebmly and annotation within Anvio
################################################################
#start screen
screen -S SO289_anvio

#load modules
module load conda

#Activate conda environment:
conda activate anvio-dev

#Define working directory
WORKDIR=/lisc/scratch/oceanography/efadeev/SO289/SO289_coassembly/ && cd $WORKDIR
PROJECT=SO289
REPO_DIR=/lisc/scratch/oceanography/efadeev/anvio_scripts

#generate default workflow
anvi-run-workflow -w metagenomics --get-default-config coassembly_config.json

#adjust the workflow file manualy and plot it
anvi-run-workflow -w metagenomics -c coassembly_config.json --save-workflow-graph

#run workflow
anvi-run-workflow -w metagenomics -c coassembly_config.json \
--additional-params --cluster-config cluster_config.json --jobs 200 \
--cluster "sbatch --output=$WORKDIR/slurm_out/%x-%j.out --partition={cluster.slurm_partition} \
--mem={cluster.mem_gb}GB --cpus-per-task={cluster.nodes} --time=1-{cluster.runtime}"

#exit screen using Ctrl+A and then D
#check status `screen -r SO289_anvio'

################################################################
#Taxonomically classify raw reads
################################################################
mkdir $WORKDIR/05_TAXONOMY

conda activate kaiju

sbatch -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT,REPO_DIR=$REPO_DIR \
--job-name "run_kaiju_fastq" $REPO_DIR/slurm_scripts/run_kaiju_fastq.sh


################################################################
#generate reference protein fasta for PD
################################################################
mkdir $WORKDIR/11_PROTEIN

sbatch -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT,REPO_DIR=$REPO_DIR \
--job-name "generate_ref_prot4PD" $REPO_DIR/slurm_scripts/run_cdhit.sh

#Import the reference databases into PD and run the MS analysis

################################################################
#Generate gene caller IDs of the idetnfied proteins 
################################################################
mkdir $WORKDIR/11_PROTEIN/R_data
#place samples_meta.txt and peptide results table into 11_PROTEIN/R_data

Rscript $REPO_DIR/slurm_scripts/detected_proteins_list.R

prot_GCIDs=$(cat 11_PROTEIN/R_data/prot_GCIDs.txt)

################################################################
#Export basic annotation of proteins from Anvio
################################################################

DATABASES=("COG20_FUNCTION" "COG20_CATEGORY")
for db in "${DATABASES[@]}"
do
sbatch -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT,db=$db \
--job-name "export_protein_annotation" $REPO_DIR/slurm_scripts/export_anvio_annotations.sh
done

"KOfam"
################################################################
#Further annotate only proteins that were identified in the MS analysis
################################################################

anvi-get-sequences-for-gene-calls --gene-caller-ids $prot_GCIDs --contigs-db $WORKDIR/03_CONTIGS/$PROJECT-contigs.db \
--get-aa-sequences --output-file $WORKDIR/11_PROTEIN/$PROJECT-detected_proteins.faa

#split the proteins fasta into chucks
mkdir $WORKDIR/11_PROTEIN/temp_fasta

awk 'BEGIN {n=0;} /^>/ {if(n%100==0){file=sprintf("./11_PROTEIN/temp_fasta/SO289-prot.faa_chunk%d",n);} print >> file; n++; next;} \
{ print >> file; }' < $WORKDIR/11_PROTEIN/$PROJECT-detected_proteins.faa

#create chunks list
ls $WORKDIR/11_PROTEIN/temp_fasta/SO289-prot.faa_chunk* > $WORKDIR/11_PROTEIN/chunk.list

#check how many chunck are there 
n_chunks=$(awk 'END { print NR }' ./11_PROTEIN/chunk.list)

################################################################
#Add taxonomy of proteins using Kaiju
################################################################
mkdir $WORKDIR/11_PROTEIN/temp_taxa

conda activate kaiju

sbatch --array=0-$n_chunks -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT,REPO_DIR=$REPO_DIR \
--job-name "run_kaiju_taxonomy" $REPO_DIR/slurm_scripts/run_kaiju_nr.sh

sbatch --array=0-$n_chunks -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT,REPO_DIR=$REPO_DIR \
--job-name "run_kaiju_taxonomy" $REPO_DIR/slurm_scripts/run_kaiju_refseq.sh

#merge chuncks

cat $WORKDIR/11_PROTEIN/temp_taxa/*_kaiju_db_refseq_2024-08-14_names.out|
awk 'BEGIN{FS=OFS="\t"}; NR==1{print "gene_callers_id","NCBI_taxID","taxa"}; \
{print $2, $3, $8}' - > $WORKDIR/11_PROTEIN/$PROJECT-detected_proteins_kaiju_taxa_refseq.tsv

cat $WORKDIR/11_PROTEIN/temp_taxa/*_kaiju_db_nr_2024-08-25_names.out|
awk 'BEGIN{FS=OFS="\t"}; NR==1{print "gene_callers_id","NCBI_taxID","taxa"}; \
{print $2, $3, $8}' - > $WORKDIR/11_PROTEIN/$PROJECT-detected_proteins_kaiju_taxa_nr.tsv

################################
#Run Interproscan (functional annotation)
################################
mkdir $WORKDIR/11_PROTEIN/temp_InterPro

#define databases 
DATABASES=("NCBIfam" "Pfam" "SignalP_GRAM_NEGATIVE" "SignalP_GRAM_POSITIVE" "PANTHER" "SUPERFAMILY" "ProSitePatterns" "TIGRFAM" "Phobius")

for db in "${DATABASES[@]}"
do
#adjust array and run sbatch
sbatch --array=0-$n_chunks -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT,db=$db --job-name "run_interproscan" $REPO_DIR/slurm_scripts/run_interproscan.sh
done

#combine chuncks of annotations
for db in "${DATABASES[@]}"
do
dbs_with_stats=("NCBIfam" "Pfam" "PANTHER" "TIGRFAM" "SUPERFAMILY")

if [[ ${dbs_with_stats[@]} =~ $db ]]; 
then 
cat $WORKDIR/11_PROTEIN/temp_InterPro/*_${db}_interpro-output.tsv |
awk 'BEGIN{FS=OFS="\t"};$9 < 0.001' -|
awk 'BEGIN{FS=OFS="\t"}; NR==1{print "gene_callers_id","source","accession","function"}; \
{print $1, $4, $5, $6}' - > $WORKDIR/11_PROTEIN/$PROJECT-detected_proteins_${db}.tsv

else 
cat $WORKDIR/11_PROTEIN/temp_InterPro/*_${db}_interpro-output.tsv |
awk 'BEGIN{FS=OFS="\t"}; NR==1{print "gene_callers_id","source","accession","function"}; \
{print $1, $4, $5, $6}' - > $WORKDIR/11_PROTEIN/$PROJECT-detected_proteins_${db}.tsv

fi
done


#get InterPro annotations
cat $WORKDIR/11_PROTEIN/temp_InterPro/*_Pfam_interpro-output.tsv |
awk 'BEGIN{FS=OFS="\t"}; NR==1{print "gene_callers_id","source","accession","function"}; \
{print $1, "InterPro", $12, $13}' - > $WORKDIR/11_PROTEIN/$PROJECT-detected_proteins_InterPro.tsv



rm -Rf $WORKDIR/11_PROTEIN/temp_InterPro
################################
#Run blastp (based on ncbi bacterial refseq annotation)
################################
mkdir $WORKDIR/11_PROTEIN/temp_blastp

#generate bacterial refseq protein database 
sbatch -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT --job-name "makeblastdb_refseq" $REPO_DIR/slurm_scripts/make_refseq_blastp_db.sh

#run blastp
sbatch --array=0-$n_chunks -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT --job-name "run_blastp" $REPO_DIR/slurm_scripts/run_blastp_refseq.sh

#merge output from all chuncks
cat $WORKDIR/11_PROTEIN/temp_blastp/*_bacterial_refseq_annotated.tsv |
awk 'BEGIN{FS=OFS="\t"};$3 > 25 && $5 < 0.001 && $6 > 50 {print $1, $2, $7, $8, $11}' - |
awk 'BEGIN{FS=OFS="\t"}; NR==1{print "gene_callers_id","blastp_acc","blastp_ann","blastp_species","blastp_taxa"}; \
{print}' - > $WORKDIR/11_PROTEIN/$PROJECT-detected_proteins_blastp_refseq.txt

################################
#Run dbcan (identify CAZymes)
################################
mkdir $WORKDIR/11_PROTEIN/temp_dbcan/

sbatch --array=0-$n_chunks -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT --job-name "run_dbcan" $REPO_DIR/slurm_scripts/run_dbcan.sh

#merge output from all chuncks
cat $WORKDIR/11_PROTEIN/temp_dbcan/*dbcan_results.tsv| 
awk 'BEGIN{OFS="\t"}; NR==1{print "gene_callers_id","EC_number","annotation"}; \
{print}' - > $WORKDIR/11_PROTEIN/$PROJECT-detected_proteins_dbcan.txt

################################
#Run MEROPS (identify peptidases)
################################
mkdir $WORKDIR/11_PROTEIN/temp_merops/

sbatch --array=0-$n_chunks -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT --job-name "run_merops" $REPO_DIR/slurm_scripts/run_merops.sh

#merge output from all chuncks
cat $WORKDIR/11_PROTEIN/temp_merops/*_merops.ann |
awk 'BEGIN{OFS="\t"}; NR==1{print "gene_callers_id","source","accession","function","e_value"}; \
{print}' - > $WORKDIR/11_PROTEIN/$PROJECT-detected_proteins_merops.txt


################################
#Run psortB (sub-cellular localization)
################################
mkdir $WORKDIR/11_PROTEIN/temp_psortb

#adjust array and run sbatch
sbatch --array=0-$n_chunks -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT --job-name "run_PSORTb_neg" $REPO_DIR/run_psortb.sh

################################
#Run deeplocpro  (sub-cellular localization)
################################
mkdir $WORKDIR/11_PROTEIN/temp_deeploc 

#adjust array and run sbatch
sbatch --array=0-$n_chunks -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT --job-name "run_deeploc" $REPO_DIR/slurm_scripts/run_deeplocpro.sh

#merge
cat $WORKDIR/11_PROTEIN/temp_deeploc/*/results_*.csv | grep -v "ACC" - | \
awk 'BEGIN{FS=","; OFS="\t"}; NR==1{print "gene_callers_id","deeploc"}; {print $2, $3}' - > $WORKDIR/11_PROTEIN/$PROJECT-detected_proteins_deeploc.txt
