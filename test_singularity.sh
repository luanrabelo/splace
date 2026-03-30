SIF_PATH=splace_v4.sif
INPUT_DIR=cp
OUTPUT_DIR=results_cp
THREADS=12

mkdir -p $OUTPUT_DIR
singularity run \
--bind ${INPUT_DIR}:/app/input_data \
--bind ${OUTPUT_DIR}:/app/output_data \
"$SIF_PATH" \
-i /app/input_data \
-o /app/output_data \
--gb-type cp \
--align \
--trimal \
--threads "$THREADS" \
--iqtree \
--feature-types CDS 