library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(GenomeInfoDb)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readxl)

# Following sc-ATAC analysis from Signac vignette: https://stuartlab.org/signac/articles/pbmc_vignette

counts <- Read10X_h5(filename = "filtered_peak_bc_matrix.h5")

metadata <- read.csv(
  file = "singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)


granges(pbmc)

peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]


library(AnnotationHub)
ah <- AnnotationHub()

# Search for the Ensembl 112 EnsDb for Homo sapiens on AnnotationHub
query(ah, "EnsDb.Hsapiens.v112")

ensdb_v112 <- ah[["AH116860"]]

annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v112)

# add the gene information to the object
Annotation(pbmc) <- annotations

# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc)

# add fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

# add blacklist ratio
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

# add reads in DNAse regions
pbmc$pct_reads_in_DNase <- pbmc$DNase_sensitive_region_fragments / pbmc$passed_filters * 100

DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 1.5, 'NS > 1.5', 'NS < 1.5')
table(pbmc$nucleosome_group)
FragmentHistogram(object = pbmc,group.by = 'nucleosome_group', region = "1-1-10000000")

VlnPlot(
  object = pbmc,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks',"pct_reads_in_DNase"),
  pt.size = 0.1,
  ncol = 3
)


#Remove outliers based on QC plots
pbmc <- subset(
  x = pbmc,
  subset = nCount_peaks > 7700 &
    nCount_peaks < 120000 &
    pct_reads_in_peaks > 40 &
    nucleosome_signal < 1.5 &
    TSS.enrichment > 4
)

#Normalization and linear dimensional reduction
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

DepthCor(pbmc)

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

DimPlot(object = pbmc, label = TRUE) + NoLegend()

#Create a gene activity matrix
gene.activities <- GeneActivity(pbmc)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'CD4'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)


#Filter barcodes based on mtDNA----
mgatk_se <- readRDS("test_scARC.rds")

# Subset to HQ cells that exist so far
mgatk_se <- mgatk_se[,colnames(pbmc)]


# Threshold based on abundance of mtDNA
pbmc$mtDNA_depth <- mgatk_se$depth

qplot(pbmc$mtDNA_depth) + scale_x_log10(breaks = c(1,10,100)) +
  geom_vline(xintercept = 10, color = "firebrick") +
  labs(x = "mean mtDNA depth/bp", y = "count") + theme_classic()

# Final filter -- mgatk, Signac, and fragments all
pbmc_filt <- subset(pbmc, mtDNA_depth >=10)

mgatk_se <- mgatk_se[,colnames(pbmc_filt)]

FilterCells(
  fragments = fragment.path,
  cells = colnames(pbmc_filt),
  outfile  = "pbmc.fragments.tsv.gz"
  )


# load mgatk output
mito.data <- ReadMGATK(dir = "mgatk_output_final/")

# create an assay
mito <- CreateAssayObject(counts = mito.data$counts)

# Subset to cell present in the scATAC-seq assay
mito <- subset(mito, cells = Cells(pbmc))

pbmc[["mito"]] <- mito

pbmc <- AddMetaData(pbmc, metadata = mito.data$depth[Cells(mito), ], col.name = "mtDNA_depth")

# Filter cells based on mtDNA depth
VlnPlot(pbmc, "mtDNA_depth", pt.size = 0.1 ) + scale_y_log10()

pbmc <- subset(pbmc, mtDNA_depth >= 10)

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 10)
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:50)
pbmc <- FindNeighbors(pbmc, reduction = "lsi", dims = 2:50)
pbmc <- FindClusters(pbmc, resolution = 0.7, algorithm = 3)

DimPlot(pbmc, label = TRUE, pt.size=2) + NoLegend()

# compute gene accessibility
DefaultAssay(pbmc) <- 'peaks'
gene.activities <- GeneActivity(pbmc)

# add to the Seurat object as a new assay
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)

pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(
  object = pbmc,
  features = c('TREM1', 'EPCAM', "PTPRC", "IL1RL1","GATA3", "KIT"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 2
)

p1 <- FeatureScatter(pbmc, "mtDNA_depth", "pct_reads_in_peaks") + ggtitle("") + scale_x_log10()
p2 <- FeatureScatter(pbmc, "mtDNA_depth", "nucleosome_signal") + ggtitle("") + scale_x_log10()

p1 + p2 + plot_layout(guides = 'collect')

variable.sites <- IdentifyVariants(pbmc, assay = "mito", refallele = mito.data$refallele)
VariantPlot(variants = variable.sites)

# Establish a filtered data frame of variants based on this processing
high.conf <- subset(
  variable.sites, subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)

high.conf[,c(1,2,5)]

pbmc <- AlleleFreq(
  object = pbmc,
  variants = high.conf$variant,
  assay = "mito"
)

pbmc[["alleles"]]

DefaultAssay(pbmc) <- "alleles"

alleles.view <- c("9520G>A", "14378T>C", "7775G>A", "204T>C")
FeaturePlot(
  object = pbmc,
  features = alleles.view,
  order = TRUE,
  cols = c("grey", "darkred"),
  pt.size = 2,
  ncol = 4
) & NoLegend()

DoHeatmap(pbmc, features = rownames(pbmc), slot = "data", disp.max = 1) +
  scale_fill_viridis_c()

#-----------------------------------------
# Read EV data from mitopore
SEV <- read_xlsx("V:/Team_LabNotebooks/TDA/Iris_Experiments/8_ATAC-Seq/8.3_2.5.2024_scATAC-Seq Trial.1/8.3_2.5.2024_EV/Analysis/Results_SEV/mitopore/mitopore_detection_threshold_0.01/Analysis/Results/disease_table.xlsx")
LEV <- read_xlsx("V:/Team_LabNotebooks/TDA/Iris_Experiments/8_ATAC-Seq/8.3_2.5.2024_scATAC-Seq Trial.1/8.3_2.5.2024_EV/Analysis/Results_LEV/mitopore/mitopore_detection_threshold_0.01/Analysis/Results/disease_table.xlsx")

LEV_server <- read.table("V:/Team_LabNotebooks/TDA/Iris_Experiments/8_ATAC-Seq/8.3_2.5.2024_scATAC-Seq Trial.1/8.3_2.5.2024_EV/Analysis/Results_LEV/mtDNA-Server2/mtDNA-Server_LEV_0.01/variants.annotated.txt", header = TRUE, sep = "\t")
SEV_server <- read.table("V:/Team_LabNotebooks/TDA/Iris_Experiments/8_ATAC-Seq/8.3_2.5.2024_scATAC-Seq Trial.1/8.3_2.5.2024_EV/Analysis/Results_SEV/mtDNA-Server2/mtDNA-Server_SEV_0.01/variants.annotated.txt", header = TRUE, sep = "\t")

SEV <- SEV %>%
  mutate(Mutation = paste0(Pos, trans))

LEV <- LEV %>%
  mutate(Mutation = paste0(Pos, trans))

# Create merged columns for the Server data
LEV_server <- LEV_server %>%
  mutate(Mutation = paste0(Pos, Ref, ">", Variant))
SEV_server <- SEV_server %>%
  mutate(Mutation = paste0(Pos, Ref, ">", Variant))

Cell <- high.conf
colnames(Cell)[3] <- c("Mutation")
Cell <- Cell[,3:6]

# Create Cell_Traced as a copy of Cell
Cell_Traced <- Cell

# Add columns for LEV and SEV matches
Cell_Traced$LEV_Match <- Cell_Traced$Mutation %in% LEV_server$Mutation
Cell_Traced$SEV_Match <- Cell_Traced$Mutation %in% SEV_server$Mutation

# Filter Cell_Traced to keep only rows with matches
Cell_Traced <- Cell_Traced[Cell_Traced$LEV_Match | Cell_Traced$SEV_Match, ]
