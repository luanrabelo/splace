SIF_PATH=splace.sif
INPUT_DIR=datasets/mitochondrial/genbank/Bufonidae
OUTPUT_DIR=results/Bufonidae
THREADS=12

mkdir -p $OUTPUT_DIR
singularity run \
--bind ${INPUT_DIR}:/app/input_data \
--bind ${OUTPUT_DIR}:/app/output_data \
"$SIF_PATH" \
-i /app/input_data \
-o /app/output_data \
--gb-type mt \
--align \
--trimal \
--threads "$THREADS" \
--iqtree