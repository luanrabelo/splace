#!/bin/bash
#SBATCH --job-name=splace
#SBATCH --output=%x_%j.out.txt
#SBATCH --error=%x_%j.err.txt
#SBATCH --partition=op1_slurm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --cpus-per-task=16
#SBATCH --time=240:00:00
#SBATCH --mail-user=luan.rabelo@pq.itv.org

printf "Begin: `date`\n";
TBEGIN=`echo "print time();" | perl`

printf "\n"

# Move to workdir
WORKDIR=$SLURM_SUBMIT_DIR

cd $WORKDIR


# Create temp dirs for proper execution
mkdir -p $WORKDIR/singularity $WORKDIR/tmp $WORKDIR/$OUTPUT_DIR

export SINGULARITY_CACHEDIR=$WORKDIR/singularity
export TMPDIR=$WORKDIR/tmp

SIF_PATH=splace.sif
INPUT_DIR=datasets/mitochondrial/genbank/Bufonidae
OUTPUT_DIR=results/Bufonidae
THREADS=12

singularity run \
--bind $WORKDIR/$INPUT_DIR:/app/input_data \
--bind $WORKDIR/$OUTPUT_DIR:/app/output_data \
"$SIF_PATH" \
-i /app/input_data \
-o /app/output_data \
--gb-type mt \
--align \
--trimal \
--threads "$THREADS" \
--iqtree

TEND=`echo "print time();" | perl`

printf "\n"
printf "End: `date`\n";
printf "Run time (s): `expr $TEND - $TBEGIN`\n";
printf "Run time (min): `expr $(( ($TEND - $TBEGIN)/60 ))`\n";