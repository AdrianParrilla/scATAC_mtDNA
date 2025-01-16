#!/bin/bash

# This script converts a bamfile to a fastq file

WORK_DIR="/mnt/smb/TDA/Iris_nextflow/1.4-1.7_Aging_ATAC-Seq/4_sinto_sample1/bams_undup"
OUT_DIR="/mnt/smb/TDA/Iris_nextflow/1.4-1.7_Aging_ATAC-Seq/6_cellranger_run2/fastqs_sample1"

for BAM in $WORK_DIR/*.bam
do
    # Extract the base name of the BAM file (without extension)
    BAM_NAME=$(basename $BAM _undup_pair.bam)
    echo "Processing $BAM_NAME..."

    bamtofastq --nthreads=10 --relaxed --reads-per-fastq=100000000 $BAM "${OUT_DIR}/$BAM_NAME"

    echo "Fastq generated for $BAM_NAME"
done