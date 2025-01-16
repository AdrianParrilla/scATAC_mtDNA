#!/bin/bash

# cellSNPs parameters
BAM=/mnt/smb/TDA/Iris_nextflow/1.4-1.7_Aging_ATAC-Seq/1_cellranger_run1/cellranger_output_sample1/outs/possorted_bam.bam
BARCODE=/mnt/smb/TDA/Iris_nextflow/1.4-1.7_Aging_ATAC-Seq/1_cellranger_run1/cellranger_output_sample1/outs/filtered_peak_bc_matrix/barcodes.tsv 
OUT_DIR_CELLSNP=/mnt/smb/TDA/Iris_nextflow/1.4-1.7_Aging_ATAC-Seq/2_cellSNP_sample1
REGION_VCF=/mnt/smb/TDA/Iris_nextflow/1.4-1.7_Aging_ATAC-Seq/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz

# The output directory of cellSNP is the input of Vireo
OUT_DIR_VIRE0=/mnt/smb/TDA/Iris_nextflow/1.4-1.7_Aging_ATAC-Seq/3_VireoSNP_sample1
n_donor=3

# The Bam file that Sinto takes is the initial possorted.bam
# donor_filter.tsv is the file generated in R
DONOR_IDS="/mnt/smb/TDA/Iris_nextflow/1.4-1.7_Aging_ATAC-Seq/3_VireoSNP_sample1/donors_filter.tsv" 
OUT_DIR="/mnt/smb/TDA/Iris_nextflow/1.4-1.7_Aging_ATAC-Seq/4_sinto_sample1"

cellsnp-lite -s $BAM -b $BARCODE -O $OUT_DIR -R $REGION_VCF -p 30 --minMAF 0.1 --minCOUNT 20 --cellTAG CB --UMItag None --gzip --genotype

vireo -c $OUT_DIR -N $n_donor -o $OUT_DIR

sinto filterbarcodes -b $BAM -c $DONOR_IDS -p 4 --outdir $OUT_DIR 