#!/bin/bash/

#First we create the masked genome to avoid mitochondrial homologous regions in the nuclear genome (step 43 of the protocol)
#Blacklist bed file was downloaded from https://github.com/caleblareau/mitoblacklist/

bedtools maskfasta -fi genome.fa -bed mtDNA_blacklist.bed -fo masked_genome.fa

# Generate a config file including this new reference, this has the format:
{
    organism: "human"
    genome: ["GRCh38"]
    input_fasta: ["/path/to/GRCh38/masked_genome.fa"]
    input_gtf: ["/path/to/gencode/Homo_sapiens.GRCh38.112.gtf"]
    non_nuclear_contigs: ["MT"]
    input_motifs: "/path/to/jaspar/motifs.pfm"
    }

#Update CellRanger-atac reference

cellranger-atac mkref --config= GRCh38_mask.config

#This will generate a new folder called GRCh38 

#Run CellRanger-atac count

cellranger-atac count --fastqs /path/to/fastqs-directory --id mgatk_out --reference /path/to/GRCh38 --localcores 60 --localmem 300 

#For mitochondrial heteroplasmy analysis, mgatk was installed in a conda environment called mito. It needs some additional dependencies.

mamba create -y -n mito openjdk r-data.table r-matrix bioconductor-genomicranges bioconductor-summarizedexperiment matplotlib

#Activate the environment and install mgatk
conda activate mito
pip install mgatk

# In order to make it work I followed suggestions from: https://github.com/caleblareau/mgatk/issues/60 and https://github.com/caleblareau/mgatk/pull/87/commits/5d9b665e00b99dc9e6b5dd382e8be592e09bee4b
# I changed source code variant_calling.py at ~/miniconda3/envs/mito/lib/python3.12/site-packages/mgatk/bin/python/variant_calling.py
# In line 159 I changed .astype(np.float) to .astype(np.float64)

#Running mgatk analysis

mgatk tenx -i /home/path/to/possorted_bam.bam -n test_scARC -o test_scARC_mgatk_V2 -b /home/path/to/filtered_peak_bc_matrix/barcodes.tsv -g 'GRCh38' -c 60 -bt CB --keep-temp-files


