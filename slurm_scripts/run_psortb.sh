#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --partition=basic
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=4:00:00

pushd $TMPDIR

git clone https://github.com/brinkmanlab/psortb_commandline_docker.git

cd psortb_commandline_docker
docker build --rm --no-cache -t psortb_commandline_docker .
sed -i "s/sudo //" psortb #remove sudo 
sed -i "s/docker run/docker run --storage-opt ignore_chown_errors=true/" psortb # resolve https://unix.stackexchange.com/a/689179
sed -i "s@brinkmanlab@docker.io/brinkmanlab@" psortb # avoid interactive choice of repo

# Exit the slurm script if a command fails
set -e

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

./psortb -i $chunk -r $TMPDIR/${file} --negative --output long

mv $TMPDIR/${file} $WORKDIR/11_PROTEIN/temp_psortb/

done

rm -rf $TMPDIR