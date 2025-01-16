#!/bin/bash

# This script removes duplicates, and unpaired reads from the bam file and sort it by name 
# Which is a requirement for bamtofastq command


# Define working directory
WORK_DIR=/mnt/smb/TDA/Iris_nextflow/1.4-1.7_Aging_ATAC-Seq/4_sinto_sample1
OUT_DIR=/mnt/smb/TDA/Iris_nextflow/1.4-1.7_Aging_ATAC-Seq/4_sinto_sample1/bams_undup

mkdir -p $OUT_DIR

# Loop through all BAM files in the directory
for BAM in $WORK_DIR/*.bam
do
    # Extract the base name of the BAM file (without extension)
    BAM_NAME=$(basename $BAM .bam)

    echo "Processing $BAM_NAME..."
    
    # Step 1: Sort BAM by read names (required for fixmate)
    samtools sort -n -o "${OUT_DIR}/${BAM_NAME}_names.bam" -O BAM $BAM

    # Step 2: fill in mate coordinates and insert size fields (required for markdup)
    samtools fixmate -m "${OUT_DIR}/${BAM_NAME}_names.bam" "${OUT_DIR}/${BAM_NAME}_fixmate.bam"

    # Step 3: Sort BAM by genomic coordinates (required for markdup)
    samtools sort -o "${OUT_DIR}/${BAM_NAME}_sorted.bam" "${OUT_DIR}/${BAM_NAME}_fixmate.bam"
    
    # Step 4: Remove duplicates
    samtools markdup -r -s "${OUT_DIR}/${BAM_NAME}_sorted.bam" "${OUT_DIR}/${BAM_NAME}_sorted_undup.bam"

    # Step 4: Remove unpaired reads/ singletons
    samtools view -h -f 2 -b "${OUT_DIR}/${BAM_NAME}_sorted_undup.bam" > "${OUT_DIR}/${BAM_NAME}_sorted_undup_pair.bam"

    # Step 5: Order again for bamtofastq
    samtools sort -n -o "${OUT_DIR}/${BAM_NAME}_undup_pair.bam" -O BAM "${OUT_DIR}/${BAM_NAME}_sorted_undup_pair.bam"

    # Clean up intermediate files to save disk space
    rm "${OUT_DIR}/${BAM_NAME}_fixmate.bam" \
        "${OUT_DIR}/${BAM_NAME}_sorted.bam" \
        "${OUT_DIR}/${BAM_NAME}_names.bam" \
        "${OUT_DIR}/${BAM_NAME}_sorted_undup.bam" \
        "${OUT_DIR}/${BAM_NAME}_sorted_undup_pair.bam"

    echo "File $BAM_NAME done"
done