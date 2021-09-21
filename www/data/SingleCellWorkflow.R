library(Seurat)
library(DropletUtils)
library(DoubletFinder)
library(clustree)
library(scater)
library(tidyverse)
library(dplyr)
library(tidyr)
library(rafalib)
library(ggplot2)
library(tvthemes)
library(ggrepel)
library(ggridges)
library(cowplot)
library(RColorBrewer)
library(viridis) 
display.brewer.all(colorblindFriendly = T)
library(ggthemes)
library(data.table)

# Set global theme and color
theme_set(theme_minimal())
# scale_colour_colorblind
#scale_colour_discrete <- scale_colour_colorblind

###########################################################################################
# 0. Input data and generate Seurat objects
###########################################################################################
#=========================================================================================
# 0.0 Set Parameters
#=========================================================================================
## Loading Cell Ranger count output
prefix <- "ProjectName"
cellRangercount <- "FilteredFeatureBcMatrix"
## Quality Control
mads <- QCMADs
raw_nCount_RNA_min <- RawCountRNAMin
raw_nCount_RNA_max <- RawCountRNAMax
raw_nFeature_RNA_min <- RawFeatureRNAMin
raw_nFeature_RNA_max <- RawFeatureRNAMax
qclog10GenesPerUMI <- Log10GenesPerUMIValue
qcmitoRatio <- MitoRatioValue
geneExpressiomCell <- GeneExpressiomAtLeastCell
removeGenesZeroCell <- "RemoveGenesThatExpressionValusEqualtoZeroInAllCell"
## Normalization and ScaleData
normalMethod <- "NormalizationMethod"
maxnfeatures <- NumberOfMostVariantGenes
## Dimensionality Reduction
clusteralgorithm <- unlist(strsplit("ClusteringAlgorithm", ","))
maxdimspca <- PCAMaxdims
perplexitystsne <- tSNEPerplexitys
mindistsumap <- UMAPMindists
## Clustering
resolutionsCluster <- c(ResolutionsForCluster)
## Marker Genes Identification
numberClusters <- ExpectedNumberOfCellClusters
topGenePlot <- TopGeneToPlotDoHeatmapFeaturePlotVlnPlotAndDotPlot
posType <- "OnlyPosMarkers"
pctmin <- MinPctMarkers
usethresh <- ThreshUseMarkers
thresholdlogfc <- LogfcThresholdMarkers
#=========================================================================================
# 0.1 Input data
#=========================================================================================
# Set root directory
dir.create(paste0(prefix, "_SingleCell_Analysis_", maxnfeatures))
setwd(paste0(prefix, "_SingleCell_Analysis_", maxnfeatures))
cat(paste0("############################# The work directory is ###########################\n", getwd(), "\n###############################################################################\n"))

# Set external sample metadata
aggregation <- read.csv("SamplesMetadata")

# Adding Gene Annotations and Mitochondrial genes
annotations <- read.delim("GenesAnnotationFile", sep = "\t", header = FALSE)
annotations <- annotations %>% separate(V1, c("GeneID1", "GeneID2", "GeneID3"), sep = "([_\\.])")
annotations <- annotations  %>%  unite(gene_name, c(GeneID1, GeneID2), sep = "-", remove = TRUE)
annotations <- annotations[ , c("gene_name", "V5")]
colnames(annotations)[2] <- 'Annotation_uniprot'

mito_genes_data <- read.csv("MitochondriaGenes")
mito_genes <- mito_genes_data$mito_genes

chlo_genes_data <- read.csv("ChloroplastGenes")
chlo_genes <- chlo_genes_data$chlo_genes

# Marker Genes
mgene <- read.delim("MarkerGenesOfKnownCellType")

# KEGG database
# Download the JSON file (this is often updated)
# https://www.genome.jp/kegg-bin/get_htext?ko00001
kegg_json <- "ko00001.json"

# Adjust the allowed object size limit in R. Example 100G
options(future.globals.maxSize = 100000 * 1024^2)

#=========================================================================================
# 0.2 Generate Seurat objects
#=========================================================================================
## Combination of multiple samples
for (file in levels(aggregation$AMLStatus)){
  seurat_data <- Read10X(data.dir = cellRangercount)
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 200, 
                                   project = file)
  assign(file, seurat_obj)
}
for (i in levels(aggregation$AMLStatus)) {
  print(get(i))
}
# Create a merged Seurat object
pscrna_seurat <- merge(as.factor(levels(aggregation$AMLStatus))[1], y = as.factor(levels(aggregation$AMLStatus))[-1], 
                       add.cell.ids = levels(aggregation$AMLStatus),
                       project = prefix)
pscrna_seurat
# Explore the metadata
head(pscrna_seurat@meta.data)
###########################################################################################
# 1. Quality control of raw scRNA-seq data
# library size and number of expressed genes
###########################################################################################
#=========================================================================================
# 1.1 Preprocessing of metadata
#=========================================================================================
# Add number of genes per UMI for each cell to metadata
pscrna_seurat$log10GenesPerUMI <- log10(pscrna_seurat$nFeature_RNA) / log10(pscrna_seurat$nCount_RNA)

# Percentage of Largest Gene (The calculation needs large memory)
# run apply over the columns (cells) and calculate what percentage of the data comes from the single most observed gene. 
# Again, having a high proportion of your data dominated by a single gene would be a concerning sign.

# pscrna_seurat$PercentLargestGene <- apply(pscrna_seurat@assays$RNA@counts, 2, function(x)(100*max(x))/sum(x))

# Compute percent mito ratio
pscrna_seurat$mitoRatio <- PercentageFeatureSet(object = pscrna_seurat, features = mito_genes)
pscrna_seurat$mitoRatio <- pscrna_seurat@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- pscrna_seurat@meta.data
# Add cell IDs and name to metadata
metadata$cells <- rownames(metadata)
metadata <- metadata %>% 
            separate(cells, c("sapmlebarcode", "order"), "-", remove = FALSE) %>% 
            separate(sapmlebarcode, c("sample", "barcode"), "_") %>% 
            unite(library_id, c(sample, order), sep = "_", remove = TRUE) %>% 
            left_join(aggregation, by = c("library_id" = "library_id")) %>% 
            select(-molecule_h5, -barcode, -AMLStatus) %>% 
            dplyr::rename(raw_orig.ident = orig.ident,
                          raw_nCount_RNA = nCount_RNA,
                          raw_nFeature_RNA = nFeature_RNA)
            
rownames(metadata) <- metadata$cells
drops <- c("cells")
metadata <- metadata[ , !(names(metadata) %in% drops)]
head(metadata)

# Add metadata back to Seurat object
pscrna_seurat@meta.data <- metadata
head(pscrna_seurat@meta.data)
# Create .RData object to load at any time
save(pscrna_seurat, file=paste0(prefix, "_seurat.RData"), compress = F)
#load(paste0(prefix, "_seurat.RData"))
#=========================================================================================
# 1.2 Evaluation of quality indicators
#=========================================================================================
# Visualize the number of cell counts per sample
pscrna_seurat$raw_orig.ident = factor(pscrna_seurat$raw_orig.ident, levels=unique(pscrna_seurat$raw_orig.ident))
pscrna_seurat$library_id = factor(pscrna_seurat$library_id, levels=unique(pscrna_seurat$library_id))

feats <- c("raw_nFeature_RNA", "raw_nCount_RNA", "log10GenesPerUMI", "mitoRatio")
pdf(paste0("1.2_", prefix, "_QC.pdf"), width = 12, height = 6)
VlnPlot(pscrna_seurat, group.by= "raw_orig.ident", features = feats, pt.size = 0, ncol = 4) + NoLegend()
VlnPlot(pscrna_seurat, group.by= "library_id", features = feats, pt.size = 0, ncol = 4) + NoLegend()

cowplot::plot_grid(ncol = 2,
                   FeatureScatter(pscrna_seurat, "raw_nCount_RNA"  , "raw_nFeature_RNA", group.by = "raw_orig.ident", pt.size = .5),
                   FeatureScatter(pscrna_seurat, "mitoRatio", "raw_nFeature_RNA", group.by = "raw_orig.ident", pt.size = .5)
)

metadata %>% 
  ggplot(aes(x=raw_orig.ident, fill=raw_orig.ident)) + 
  geom_bar() +
  theme_classic() +
  scale_fill_brewer(palette="Dark2") +
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle(paste0("Number of cell counts per sample"))

metadata %>% 
  ggplot(aes(x=library_id, fill=library_id)) + 
  geom_bar() +
  theme_classic() +
  scale_fill_brewer(palette="Dark2") +
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle(paste0("Number of cell counts per sample replicates"))

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=raw_orig.ident, x=raw_nCount_RNA, fill= raw_orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("log 10 Cell density") +
  scale_fill_brewer(palette="Dark2") +
  geom_vline(xintercept = 500) +
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle(paste0("Number UMIs/transcripts per cell"))

# Visualize Genes detected per cell
metadata %>% 
  ggplot(aes(color=raw_orig.ident, x=raw_nFeature_RNA, fill= raw_orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("log 10 Cell density") +
  scale_fill_brewer(palette="Dark2") +
  geom_vline(xintercept = 500) +
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle(paste0("Number Gene per cell"))

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=raw_nCount_RNA, y=raw_nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~ raw_orig.ident) + 
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle(paste0("UMIs vs genes"))

metadata %>%
  arrange(raw_orig.ident) %>%
  ggplot(aes(raw_nCount_RNA,raw_nFeature_RNA,colour=raw_orig.ident)) + 
  geom_point() + 
  ggtitle(paste0("UMIs vs genes")) +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  theme_classic()

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=raw_orig.ident, x=mitoRatio, fill=raw_orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2) +
  ggtitle(paste0("Mitochondrial counts ratio"))

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x = log10GenesPerUMI, color = raw_orig.ident, fill = raw_orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

#FeatureScatter(pscrna_seurat, feature1 = "nCount_RNA", feature2 = "PercentLargestGene", group.by = "orig.ident")
dev.off()
#=========================================================================================
# 1.3 MADs(scater) filter
#=========================================================================================
library(DropletUtils)
library(scater)
# The ident column is occupied when the SCE object is generated
pscrna_seurat.sce <- as.SingleCellExperiment(pscrna_seurat)
pscrna_seurat.sce <- calculateQCMetrics(pscrna_seurat.sce)
colnames(colData(pscrna_seurat.sce))
colnames(rowData(pscrna_seurat.sce))
libsize.drop <- isOutlier(pscrna_seurat.sce$total_counts, nmads=mads, type="both", log=TRUE)
feature.drop <- isOutlier(pscrna_seurat.sce$total_features_by_counts, nmads=mads, type="both", log=TRUE)
pscrna_seurat.sce <- pscrna_seurat.sce[,!(libsize.drop | feature.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), Remaining=ncol(pscrna_seurat.sce))

pdf(paste0("1.3_", prefix, "_scater_QC.pdf"), width = 12, height = 6)
plot_grid(plotColData(pscrna_seurat.sce, y = "raw_nFeature_RNA",x = "raw_orig.ident",colour_by = "raw_orig.ident"),
          plotColData(pscrna_seurat.sce, y = "raw_nCount_RNA", x = "raw_orig.ident",colour_by = "raw_orig.ident"), ncol = 2)

plotExpression(pscrna_seurat.sce, features = filter(mgene, name=="metaphloem sieve element")$gene, x = "raw_orig.ident", colour_by= "raw_orig.ident") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Compute the relative expression of each gene per cell
rel_expression <- t( t(counts(pscrna_seurat.sce)) / Matrix::colSums(counts(pscrna_seurat.sce))) * 100
most_expressed <- sort(Matrix::rowSums( rel_expression ),T)[20:1] / ncol(pscrna_seurat.sce)
par(mfrow=c(1,1),mar=c(4,8,1,1))
boxplot( as.matrix(t(rel_expression[names(most_expressed),])),cex=.1, las=1, xlab="% total count per cell",col=scales::hue_pal()(20)[20:1],horizontal=TRUE)


par(mfrow=c(2,3), mar=c(5, 4, 1, 1), bty="n")
hist(log10(pscrna_seurat.sce$total_counts), xlab="log10(Library sizes)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(log10(pscrna_seurat.sce$total_features_by_counts), xlab="log10(# of expressed genes)", 
     main="", breaks=20, col="grey80", ylab="Number of cells")
smoothScatter(log10(pscrna_seurat.sce$total_counts), log10(pscrna_seurat.sce$total_features_by_counts), 
              xlab="log10(Library sizes)", ylab="log10(# of expressed genes)", 
              nrpoints=500, cex=0.5)
hist(log10(rowData(pscrna_seurat.sce)$mean_counts+1e-6), col="grey80",  main="", 
     breaks=40, xlab="log10(ave # of UMI + 1e-6)")
hist(log10(rowData(pscrna_seurat.sce)$n_cells_by_counts+1), col="grey80", main="", 
     breaks=40, xlab="log10(# of expressed cells + 1)")
smoothScatter(log10(rowData(pscrna_seurat.sce)$mean_counts+1e-6), 
              log10(rowData(pscrna_seurat.sce)$n_cells_by_counts + 1), 
              xlab="log10(ave # of UMI + 1e-6)", 
              ylab="log10(# of expressed cells + 1)")
dev.off()

filtered_pscrna_seurat <- as.Seurat(pscrna_seurat.sce)
head(filtered_pscrna_seurat@meta.data)
#=========================================================================================
# 1.4 Seurat filter
#=========================================================================================
# Specific genes were filtered out according to 1.3 boxplot
counts <- GetAssayData(filtered_pscrna_seurat, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('genein1', 'genein2', 'genein3'))),]
filtered_pscrna_seurat <- subset(filtered_pscrna_seurat, features = rownames(counts))
filtered_pscrna_seurat <- subset(x = filtered_pscrna_seurat, 
                                 subset= (nCount_RNA >= raw_nCount_RNA_min) &
                                         (nCount_RNA <= raw_nCount_RNA_max) &
                                         (nFeature_RNA >= raw_nFeature_RNA_min) &
                                         (nFeature_RNA <= raw_nFeature_RNA_max) &
                                         (log10GenesPerUMI > qclog10GenesPerUMI) &
                                         (mitoRatio < qcmitoRatio)
                                )

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_pscrna_seurat, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than geneExpressiomCell TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= geneExpressiomCell
# Only keeping those genes expressed in more than geneExpressiomCell cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_pscrna_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_pscrna_seurat@meta.data)
#=========================================================================================
# 1.5 Detection and filtration of diploids
#=========================================================================================
library(DoubletFinder)
doubletfinder_seurat <- filtered_pscrna_seurat
doubletfinder_seurat = NormalizeData(doubletfinder_seurat)
doubletfinder_seurat <- FindVariableFeatures(doubletfinder_seurat, verbose = T)
doubletfinder_seurat <- ScaleData(doubletfinder_seurat, vars.to.regress = c("nFeature_RNA", "mitoRatio"), selection.method = "vst", nfeatures = maxnfeatures, verbose = T)
doubletfinder_seurat <- RunPCA(doubletfinder_seurat, npcs = 50, verbose = T)
doubletfinder_seurat <- FindNeighbors(doubletfinder_seurat, dims = 1:20)
doubletfinder_seurat <- FindClusters(doubletfinder_seurat, resolution = 0.5)
doubletfinder_seurat <- RunUMAP(doubletfinder_seurat, dims = 1:30, verbose = T)

## pK Identification
sweep.res <- paramSweep_v3(doubletfinder_seurat, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
save(bcmvn, file=paste0(prefix, "_doubletfinder_pK_Identification.RData"), compress = F)
#load(paste0(prefix, "_seurat.RData"))
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
pdf(paste0("1.4_", prefix, "_doubletfinder_optimize_parameters.pdf"), width = 10, height = 6)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las = 2)
dev.off()

# define the expected number of doublet cellscells
nExp <- round(ncol(doubletfinder_seurat) * 0.04)  # expect 4% doublets

doubletfinder <- doubletFinder_v3(doubletfinder_seurat, pN = 0.25, pK = pk_good, nExp = nExp, PCs = 1:30, sct = FALSE)
# name of the DF prediction can change, so extract the correct column name.
DF.name <- colnames(doubletfinder@meta.data)[grepl("DF.classification", colnames(doubletfinder@meta.data))]
pdf(paste0("1.4_", prefix, "_doubletfinder.pdf"), width = 10, height = 6)
cowplot::plot_grid(ncol = 2, 
                   DimPlot(doubletfinder, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(doubletfinder, group.by = DF.name) + NoAxes()
                   )
pv <- VlnPlot(doubletfinder, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
print(pv)
dev.off()

singlet <- pv$data %>% filter(ident=="Singlet")
doublet <- pv$data %>% filter(ident=="Doublet")

if (round(mean(doublet$nFeature_RNA) / mean(singlet$nFeature_RNA)) >= 2) {
  print("去除双胞体")
  filtered_pscrna_seurat <- doubletfinder[, doubletfinder@meta.data[, DF.name] == "Singlet"]
} else {
  print("未检测到大量双胞体的存在")
}
dim(filtered_pscrna_seurat)
head(filtered_pscrna_seurat@meta.data)
#=========================================================================================
# 1.6 Processing meta.data
#=========================================================================================
filtered_metadata <- filtered_pscrna_seurat@meta.data
filtered_drops <- c("ident")
filtered_metadata <- filtered_metadata[ , !(names(filtered_metadata) %in% filtered_drops)]
filtered_pscrna_seurat@meta.data <- filtered_metadata
head(filtered_pscrna_seurat@meta.data)
filtered_pscrna_seurat$orig.ident = factor(filtered_pscrna_seurat$orig.ident, levels=unique(filtered_pscrna_seurat$orig.ident))
filtered_pscrna_seurat$library_id = factor(filtered_pscrna_seurat$library_id, levels=unique(filtered_pscrna_seurat$library_id))
filtered_pscrna_seurat
#=========================================================================================
# 1.7 Reassess quality indicators
#=========================================================================================
filtered_pscrna_seurat.sce <- as.SingleCellExperiment(filtered_pscrna_seurat)

filtered_metadata <- filtered_pscrna_seurat@meta.data
filtered_feats <- c("nFeature_RNA", "nCount_RNA", "log10GenesPerUMI", "mitoRatio")
# Visualize the number of cell counts per sample
pdf(paste0("1.7_", prefix, "_filtered_QC.pdf"), width = 12, height = 6)
VlnPlot(filtered_pscrna_seurat, group.by= "orig.ident", features = filtered_feats, pt.size = 0, ncol = 4) + NoLegend()
VlnPlot(filtered_pscrna_seurat, group.by= "library_id", features = filtered_feats, pt.size = 0, ncol = 4) + NoLegend()

cowplot::plot_grid(ncol = 2,
                   FeatureScatter(filtered_pscrna_seurat, "nCount_RNA"  , "nFeature_RNA", group.by = "orig.ident", pt.size = .5),
                   FeatureScatter(filtered_pscrna_seurat, "mitoRatio", "nFeature_RNA", group.by = "orig.ident", pt.size = .5)
)
filtered_metadata %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  scale_fill_brewer(palette="Dark2") +
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle(paste0("Number of cell counts per sample (Filtered)"))

filtered_metadata %>% 
  ggplot(aes(x=library_id, fill=library_id)) + 
  geom_bar() +
  theme_classic() +
  scale_fill_brewer(palette="Dark2") +
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle(paste0("Number of cell counts per sample replicates"))

# Visualize the number UMIs/transcripts per cell
filtered_metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("log 10 Cell density") +
  scale_fill_brewer(palette="Dark2") +
  geom_vline(xintercept = 500) +
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle(paste0("Number UMIs/transcripts per cell (Filtered)"))

# Visualize Genes detected per cell
filtered_metadata %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("log 10 Cell density") +
  scale_fill_brewer(palette="Dark2") +
  geom_vline(xintercept = 500) +
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle(paste0("Number Gene per cell (Filtered)"))

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
filtered_metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~ orig.ident) + 
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle(paste0("UMIs vs genes (Filtered)"))

filtered_metadata %>%
  arrange(orig.ident) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=orig.ident)) + 
  geom_point() + 
  ggtitle(paste0("UMIs vs genes (Filtered)")) +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  theme_classic()

# Visualize the distribution of mitochondrial gene expression detected per cell
filtered_metadata %>% 
  ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2) +
  ggtitle(paste0("Mitochondrial counts ratio (Filtered)"))

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
filtered_metadata %>%
  ggplot(aes(x = log10GenesPerUMI, color = orig.ident, fill = orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

#Compute the relative expression of each gene per cell
rel_expression <- t( t(counts(filtered_pscrna_seurat.sce)) / Matrix::colSums(counts(filtered_pscrna_seurat.sce))) * 100
most_expressed <- sort(Matrix::rowSums( rel_expression ),T)[50:1] / ncol(filtered_pscrna_seurat.sce)
par(mfrow=c(1,1),mar=c(4,8,1,1))
boxplot(as.matrix(t(rel_expression[names(most_expressed),])),cex=.1, las=1, xlab=paste0("% total count per cell (Filtered)"), col=scales::hue_pal()(50)[50:1],horizontal=TRUE)

dev.off()
#=========================================================================================
# 1.8 Preserve the filtered cells
#=========================================================================================
# Create .RData object to load at any time
save(filtered_pscrna_seurat, file=paste0(prefix, "_filtered_seurat.RData"), compress = F)
#load(paste0(prefix, "_filtered_seurat.RData"))
###########################################################################################
# 2. Normalization and scale
###########################################################################################
norm_pscrna_seurat <- filtered_pscrna_seurat
norm_pscrna_seurat
#=========================================================================================
# 2.1 Raw distributions of their expression values
#=========================================================================================
# Split seurat object by sample (orig.ident)
split_pscrna_seurat_raw <- SplitObject(norm_pscrna_seurat, split.by = "orig.ident")
split_pscrna_seurat_raw <- split_pscrna_seurat_raw[unique(norm_pscrna_seurat$orig.ident)]

# raw distributions of their expression values
p_raw <- list()
pc_raw <- list()
for (i in unique(norm_pscrna_seurat$orig.ident)) {
  p_raw[[i]] <- as.tibble(
    split_pscrna_seurat_raw[[i]]@assays$RNA@data[,1:100]
  ) %>%
    pivot_longer(
      cols=everything(),
      names_to="cell",
      values_to="expression"
    ) %>%
    ggplot(aes(x=expression, group=cell)) +
    geom_density() +
    theme_classic() +
    coord_cartesian(ylim=c(0,1), xlim=c(0,3)) +
    ggtitle(paste0(i))
  ##########################################################
  pc_raw[[i]] <- tibble(
    pc95 = apply(split_pscrna_seurat_raw[[i]][["RNA"]]@data,2,quantile,0.95),
    measured = apply(split_pscrna_seurat_raw[[i]][["RNA"]]@data,2,function(x)(100*sum(x!=0))/length(x))
  ) %>% 
    ggplot(aes(x=measured,y=pc95))+
    geom_point()+
    theme_classic() +
    ggtitle(paste0("raw of ",i))
}
pdf(paste0("2.1_", prefix, "_Distributions_Expression_Values_raw.pdf"), width = 12, height = 6)
CombinePlots(plots = p_raw, nrow=2, legend="none")
CombinePlots(plots = pc_raw, nrow=2, legend="none")
dev.off()
#=========================================================================================
# 2.2 Normalization
#=========================================================================================
if("LogNormalize" %in% normalMethod) {
  print("Select LogNormalize to normalization")
  norm_pscrna_seurat_log <- NormalizeData(norm_pscrna_seurat, normalization.method = "LogNormalize", verbose = TRUE)
  
  # distributions of their expression values
  split_pscrna_seurat_log <- SplitObject(norm_pscrna_seurat_log, split.by = "orig.ident")
  split_pscrna_seurat_log <- split_pscrna_seurat_log[unique(norm_pscrna_seurat_log$orig.ident)]
  
  p_log <- list()
  pc_log <- list()
  for (i in unique(norm_pscrna_seurat_log$orig.ident)) {
    p_log[[i]] <- as.tibble(
      split_pscrna_seurat_log[[i]]@assays$RNA@data[,1:100]
    ) %>%
      pivot_longer(
        cols=everything(),
        names_to="cell",
        values_to="expression"
      ) %>%
      ggplot(aes(x=expression, group=cell)) +
      geom_density() +
      theme_classic() +
      coord_cartesian(ylim=c(0,1), xlim=c(0,3)) +
      ggtitle(paste0(i))
    ##########################################################
    pc_log[[i]] <- tibble(
      pc95 = apply(split_pscrna_seurat_log[[i]][["RNA"]]@data,2,quantile,0.95),
      measured = apply(split_pscrna_seurat_log[[i]][["RNA"]]@data,2,function(x)(100*sum(x!=0))/length(x))
    ) %>% 
      ggplot(aes(x=measured,y=pc95))+
      geom_point()+
      theme_classic() +
      ggtitle(paste0("Normalisation of ",i))
  }
  pdf(paste0("2.2_", prefix, "_Distributions_Expression_Values_LogNormalize.pdf"), width = 15)
  CombinePlots(plots = p_log, nrow=2, legend="none")
  CombinePlots(plots = pc_log, nrow=2, legend="none")
  dev.off()
  norm_pscrna_seurat <- norm_pscrna_seurat_log
} else if ("CLR" %in% normalMethod) {
  print("Select CLR to normalization")
  norm_pscrna_seurat_clr <- NormalizeData(norm_pscrna_seurat, normalization.method = "CLR", verbose = TRUE)
  
  # distributions of their expression values
  split_pscrna_seurat_clr <- SplitObject(norm_pscrna_seurat_clr, split.by = "orig.ident")
  split_pscrna_seurat_clr <- split_pscrna_seurat_clr[unique(norm_pscrna_seurat_clr$orig.ident)]
  
  p_clr <- list()
  pc_clr <- list()
  for (i in unique(norm_pscrna_seurat_clr$orig.ident)) {
    p_clr[[i]] <- as.tibble(
      split_pscrna_seurat_clr[[i]]@assays$RNA@data[,1:100]
    ) %>%
      pivot_longer(
        cols=everything(),
        names_to="cell",
        values_to="expression"
      ) %>%
      ggplot(aes(x=expression, group=cell)) +
      geom_density() +
      theme_classic() +
      coord_cartesian(ylim=c(0,1), xlim=c(0,3)) +
      ggtitle(paste0(i))
    ##########################################################
    pc_clr[[i]] <- tibble(
      pc95 = apply(split_pscrna_seurat_clr[[i]][["RNA"]]@data,2,quantile,0.95),
      measured = apply(split_pscrna_seurat_clr[[i]][["RNA"]]@data,2,function(x)(100*sum(x!=0))/length(x))
    ) %>% 
      ggplot(aes(x=measured,y=pc95))+
      geom_point()+
      theme_classic() +
      ggtitle(paste0("Normalisation of ",i))
  }
  pdf(paste0("2.2_", prefix, "_Distributions_Expression_Values_CLR.pdf"), width = 15)
  CombinePlots(plots = p_clr, nrow=2, legend="none")
  CombinePlots(plots = pc_clr, nrow=2, legend="none")
  dev.off()
  norm_pscrna_seurat <- norm_pscrna_seurat_clr
} else {
  print("At lest select one method")
}
#=========================================================================================
# 2.3 Identification of highly variable genes
#=========================================================================================
norm_pscrna_seurat <- FindVariableFeatures(norm_pscrna_seurat, selection.method = "vst", nfeatures = maxnfeatures)

pdf(paste0("2.3_", prefix, "_identied_variable_genes.pdf"), width = 10)
top10 <- head(VariableFeatures(norm_pscrna_seurat), 10)
vfp1 <- VariableFeaturePlot(norm_pscrna_seurat)
vfp1 <- LabelPoints(plot = vfp1, points = top10, repel = TRUE)
vfp1
dev.off()
#=========================================================================================
# 2.4 Scale
#=========================================================================================
norm_pscrna_seurat <- ScaleData(norm_pscrna_seurat, vars.to.regress = c("percent_mito", "nFeature_RNA"))
#=========================================================================================
# 2.5 Save integrated seurat object
#=========================================================================================
saveRDS(norm_pscrna_seurat, paste0(prefix, "_normalized_seurat.rds"), compress = F)
#readRDS(paste0(prefix, "_normalized_seurat.rds"))
save(norm_pscrna_seurat, file=paste0(prefix, "_normalized_seurat.RData"), compress = F)
#load(paste0(prefix, "_normalized_seurat.RData"))
###########################################################################################
# 3. PCA dimensionality reduction analysis was used to determine the maximum number of PCs for subsequent clustering analysis
# Some sources of technical differences are added to the higher PC, so the choice of PC is more important
# The more PCs you choose, the more differences you consider when performing clustering, but it takes more time to perform clustering
###########################################################################################
# Run PCA, Default calculation of 50 pcs
norm_pscrna_seurat <- RunPCA(object = norm_pscrna_seurat, features=VariableFeatures(norm_pscrna_seurat))
# Plot PCA
# Look at the top 40 PC components
p <- list()
for (i in 1:40) {
  p[[i]] <- VizDimLoadings(norm_pscrna_seurat, dims = i, ncol = 1) +
    theme_minimal(base_size = 8) +
    ggtitle(paste0("PCA loadings of PC ",i))
}

pdf(paste0("3_", prefix, "_PCA.pdf"), width = 30, height = 30)
DimPlot(norm_pscrna_seurat, reduction="pca", split.by = "orig.ident", group.by = "orig.ident")
DimPlot(norm_pscrna_seurat, reduction="pca", split.by = "library_id", group.by = "orig.ident")
DimPlot(norm_pscrna_seurat, reduction="pca", group.by = "orig.ident")
CombinePlots(plots = p, nrow=4)
dev.off()

# To identify the available dimensions of the dataset, jackdraw operation requires large memory
# The dimensions above the dotted line are available dimensions. You can also adjust the dims parameters to draw all PCA views
norm_pscrna_seurat_JS <- JackStraw(object = norm_pscrna_seurat, num.replicate = 100, dims = 40)
norm_pscrna_seurat_JS <- ScoreJackStraw(object = norm_pscrna_seurat_JS, dims = 1:40)
pdf(paste0("3_", prefix, "_JackStraw.pdf"), width = 20, height = 10)
JackStrawPlot(object = norm_pscrna_seurat_JS, dims = 1:40)
dev.off()

# Elbow plot explores PCs to see which ones have obvious differences in structure
pdf(paste0("3_", prefix, "_ElbowPlot.pdf"), width = 15)
ElbowPlot(norm_pscrna_seurat, ndims = 50)
dev.off()

# Calculate the maximum dims value
# The principal component only contributes 5% of the standard deviation, while the cumulative contribution of the principal component is 90% of the standard deviation
pct <- norm_pscrna_seurat[["pca"]]@stdev / sum(norm_pscrna_seurat[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
# Evaluate points where the percentage change between consecutive PCs is less than 0.1%
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
cat(paste0("The minimum of the two calculation is: ", pcs, "\n"))

pdf(paste0("3_", prefix, "_ElbowPlot2.pdf"), width = 15)
par(mar=c(1,1,1,1))
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

# Visualization of highly variable genes in selected PC by thermography
pdf(paste0("3_", prefix, "_DimHeatmap.pdf"), width=200, height=100)
par(mar=c(1,1,1,1))
DimHeatmap(norm_pscrna_seurat, dims=1:50, cells=500, balanced = TRUE)
dev.off()

# Observe the results of jackstrawplot, elbowplot and dimheatmap to determine the maximum number of PCs for subsequent clustering analysis [modified]
# maxdims <- pcs
maxdims <- maxdimspca
# Printing out the most variable genes driving PCs
sink(paste0("3.1_", prefix, "_most_variable_genes.txt"))
print(x = norm_pscrna_seurat[["pca"]],
      dims = 1:maxdims,
      nfeatures = 5)
sink()
###########################################################################################
# 4. Clustering and Visualizing
###########################################################################################
#=========================================================================================
# 4.1 Clustering
#=========================================================================================
cluster_pscrna_seurat <- norm_pscrna_seurat
cat(paste0("################# Begin Clustering ... #################\n"))
cluster_pscrna_seurat <- FindNeighbors(cluster_pscrna_seurat, dims=1:maxdims)
# Larger resolution values give larger clusters, smaller resolution values gives smaller clusters
resolutions <- resolutionsCluster
cluster_pscrna_seurat <- FindClusters(cluster_pscrna_seurat, resolution = resolutions)

# investigate how many clusters each resolution produces
cat("################# Stats how many clusters each resolution produces #################\n")
sapply(grep("res",colnames(cluster_pscrna_seurat@meta.data),value = TRUE),
       function(x) length(unique(cluster_pscrna_seurat@meta.data[,x])))
#=========================================================================================
# 4.2 Visualizing cell clusters
#=========================================================================================
# Explore resolutions
cat("################# Explore resolutions... #################\n")
head(cluster_pscrna_seurat@meta.data)
# clusters are stored in the “seurat_clusters” metadata
library(clustree)
clus.tree.out <- clustree(cluster_pscrna_seurat) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
pdf(paste0("4.2_", prefix, "_clustree_resolution.pdf"), width = 20)
print(clus.tree.out)
dev.off()

if("PCA" %in% clusteralgorithm) {
  print("Select PCA to Clustering and Visualizing")
  #-----------------------------------------------------------------------------------------
  # 4.2.1. PCA
  #-----------------------------------------------------------------------------------------
  pdf(paste0("4.2.1_", prefix, "_PCA_Clustering.pdf"), width = 20, height = 20)
  DimPlot(cluster_pscrna_seurat, reduction="pca", label = TRUE) + ggtitle("PC1 vs PC2 with Clusters")
  DimPlot(cluster_pscrna_seurat, reduction="pca", dims=c(4,9), label=TRUE) + ggtitle("PC4 vs PC9 with Clusters")
  dev.off()
} else if ("tSNE" %in% clusteralgorithm) {
  #-----------------------------------------------------------------------------------------
  # 4.2.2. tSNE
  #-----------------------------------------------------------------------------------------
  saved.seed <- 8888
  set.seed(saved.seed)
  perplexitys <- perplexitystsne
  resolutionsver <- paste0('cluster_pscrna_seurat_tsen_', resolutionsCluster)
  
  cluster_pscrna_seurat_tsen <- list()
  cluster_pscrna_seurat_tsen_name <- as.list(outer(resolutionsver, perplexitys, paste, sep="_"))
  pdf(paste0("4.2.2_", prefix, "_tSNE_perplexity.pdf"), width = 12, height = 10)
  for (i in 1:length(cluster_pscrna_seurat_tsen_name)) {
    resolution <- as.character(as.list(strsplit(as.character(cluster_pscrna_seurat_tsen_name[i]), split='_')[[1]])[5])
    myperplexity <- as.numeric(as.character(as.list(strsplit(as.character(cluster_pscrna_seurat_tsen_name[i]), split='_')[[1]])[6]))
    Idents(object = cluster_pscrna_seurat) <- paste0("RNA_snn_res.", resolution)
    cat(paste0("############## Stats of cells in each cluster of each sample at resolutions ", resolution, " ##############\n"))
    cat(paste0("############## At this resolution, each cluster corresponds to the number of cells in each sample ##############\n"))
    print(table(Idents(cluster_pscrna_seurat), cluster_pscrna_seurat$orig.ident))
    cat(paste0("############## Set to resolutions of ", resolution, " and perplexity ", myperplexity, " ##############\n"))
    cluster_pscrna_seurat_tsen[[i]] <- RunTSNE(cluster_pscrna_seurat,
                                               dims=1:maxdims,
                                               seed.use = saved.seed, 
                                               perplexity = myperplexity
    )
    print(cluster_pscrna_seurat_tsen[[i]])
    ## get tsen data
    phe <- data.frame(cell=rownames(cluster_pscrna_seurat_tsen[[i]]@meta.data), cluster =cluster_pscrna_seurat_tsen[[i]]@meta.data$seurat_clusters)
    tsne_pos <- Embeddings(cluster_pscrna_seurat_tsen[[i]], 'tsne')
    dat <- cbind(tsne_pos, phe)
    cat(paste0("############## Check tsen data with resolutions of ", resolution, " and perplexity ", myperplexity, " ##############\n"))
    print(head(dat))
    ## plot
    cat(paste0("############## Plot tsen with resolutions of ", resolution, " and perplexity ", myperplexity, " ##############\n"))
    p <- DimPlot(cluster_pscrna_seurat_tsen[[i]], reduction = "tsne", pt.size = 1, label = TRUE, label.size = 7) + ggtitle(paste0("tSNE with Perplexity ", myperplexity, " resolutions ", resolution, prefix))
    p2 <- DimPlot(cluster_pscrna_seurat_tsen[[i]], reduction = "tsne", pt.size = 1, label = TRUE, label.size = 7, split.by ='orig.ident') + ggtitle(paste0("tSNE with Perplexity ", myperplexity, " resolutions ", resolution, prefix))
    p3 <- DimPlot(cluster_pscrna_seurat_tsen[[i]], reduction = "tsne", pt.size = 1, label = TRUE, label.size = 7, split.by ='library_id') + ggtitle(paste0("tSNE with Perplexity ", myperplexity, " resolutions ", resolution, prefix))
    print(p)
    print(p2)
    print(p3)
  }
  dev.off()
} else if ("UMAP" %in% clusteralgorithm) {
  print("Select UMAP to Clustering and Visualizing")
  #-----------------------------------------------------------------------------------------
  # 4.2.3. UMAP可视化
  #-----------------------------------------------------------------------------------------
  mindists <- mindistsumap
  preresolutions <- paste0('cluster_pscrna_seurat_umap_', resolutionsCluster)
  
  cluster_pscrna_seurat_umap <- list()
  cluster_pscrna_seurat_umap_name <- as.list(outer(preresolutions, mindists, paste, sep="_"))
  pdf(paste0("4.2.3_", prefix, "_UMAP_Clustering.pdf"), width = 15, height = 10)
  for (i in 1:length(cluster_pscrna_seurat_umap_name)) {
    myresolution <- as.character(as.list(strsplit(as.character(cluster_pscrna_seurat_umap_name[i]), split='_')[[1]])[5])
    mymindist <- as.numeric(as.character(as.list(strsplit(as.character(cluster_pscrna_seurat_umap_name[i]), split='_')[[1]])[6]))
    Idents(object = cluster_pscrna_seurat) <- paste0("RNA_snn_res.", myresolution)
    cat(paste0("############## Set to resolutions of ", myresolution, " min.dist ", mymindist, " ##############\n"))
    print(table(Idents(cluster_pscrna_seurat), cluster_pscrna_seurat$orig.ident))
    cluster_pscrna_seurat_umap[[i]] <- RunUMAP(cluster_pscrna_seurat,
                                               dims = 1:maxdims,
                                               min.dist = mymindist,
                                               reduction = "pca")
    cat(paste0("############## Stat resolutions of ", myresolution, " min.dist ", mymindist, " ##############\n"))
    pheumap <- data.frame(cell = rownames(cluster_pscrna_seurat_umap[[i]]@meta.data),cluster = cluster_pscrna_seurat_umap[[i]]@meta.data$seurat_clusters)
    umap_pos <- Embeddings(cluster_pscrna_seurat_umap[[i]], 'umap')
    datumap <- cbind(umap_pos,pheumap)
    print(head(datumap))
    p1 <- DimPlot(cluster_pscrna_seurat_umap[[i]], pt.size=0.5, reduction = "umap", label = TRUE, repel = TRUE, label.size = 6) + 
      ggtitle(paste0("UMAP with resolution ", myresolution, " min.dist ", mymindist))
    
    p2 <- DimPlot(cluster_pscrna_seurat_umap[[i]], pt.size=0.5, reduction = "umap", label = TRUE, repel = TRUE, split.by = "orig.ident") #ggrepel::geom_text_repel(data = mgene, aes(label = genes)) + ggtitle(paste0("UMAP with resolution ", myresolution))
    
    p3 <- DimPlot(cluster_pscrna_seurat_umap[[i]], pt.size=0.5, reduction = "umap", label = TRUE, repel = TRUE, split.by = "library_id") + ggtitle(paste0("UMAP with resolution ", myresolution))
    
    print(p1)
    print(p2)
    print(p3)
  }
  dev.off()
} else {
  print("At lest select one")
}
#=========================================================================================
# 4.3 Select the appropriate resolution
#=========================================================================================
# Select the appropriate resolution based on Expected Number Of Cell Clusters and Clusters Tree in 4.2 section
# And assign the final classification to the ident column [modify]
# Here, we select the resolution 0.4
cluster_pscrna_seurat@meta.data$ident <- cluster_pscrna_seurat@meta.data$RNA_snn_res.0.4
#=========================================================================================
# 4.4 Clustering quality control
#=========================================================================================
#-----------------------------------------------------------------------------------------
# 4.4.1. Clustering by sample
#-----------------------------------------------------------------------------------------
# UMAP of cells in each cluster by orig.ident
pdf(paste0("4.4.1_", prefix, "_UMAP_Clustering_QC.pdf"), width = 30, height = 10)
DimPlot(cluster_pscrna_seurat,
        label = TRUE,
        pt.size = 1,
        split.by = "orig.ident")  + NoLegend()
dev.off()
# UMAP of cells in each cluster by condition
pdf(paste0("4.4.1_", prefix, "_UMAP_Clustering_QC2.pdf"), width = 20, height = 10)
DimPlot(cluster_pscrna_seurat,
        label = TRUE,
        pt.size = 1,
        split.by = "library_id")  + NoLegend()
dev.off()
#-----------------------------------------------------------------------------------------
# 4.4.2. Exploring the consistency of marker genes in different clusters
#-----------------------------------------------------------------------------------------
for (i in 3:length(unique(mgene$name))) {
  #for (i in c("metaphloem sieve element")) {
  cat(paste0("############### ", unique(mgene$name)[i], " ###############\n"))
  try(FeaturePlot(cluster_pscrna_seurat, reduction = "umap", features = filter(mgene, name == unique(mgene$name)[i])$gene, ncol = 3, order = T, label = TRUE, label.size = 7) + labs(subtitle = paste0("UMAP with maker ", unique(mgene$name)[i])))
  ggsave(paste0("4.4.2_", prefix, "_", unique(mgene$name)[i], "_markers_UMAP.pdf"), width = 30, height = 30, limitsize = FALSE)
  try(FeaturePlot(cluster_pscrna_seurat, reduction = "tsne", features = filter(mgene, name == unique(mgene$name)[i])$gene, ncol = 3, order = T, label = TRUE, label.size = 7) + labs(subtitle = paste0("tSNE with maker ", unique(mgene$name)[i])))
  ggsave(paste0("4.4.2_", prefix, "_", unique(mgene$name)[i], "_markers_tSNE.pdf"), width = 30, height = 30, limitsize = FALSE)
}
#-----------------------------------------------------------------------------------------
# 4.4.3. Finally, PCA, tsne and umap visualization were determined
#-----------------------------------------------------------------------------------------
cluster_pscrna_seurat$orig.ident = factor(cluster_pscrna_seurat$orig.ident, levels=unique(cluster_pscrna_seurat$orig.ident))
cluster_pscrna_seurat$library_id = factor(cluster_pscrna_seurat$library_id, levels=unique(cluster_pscrna_seurat$library_id))
pdf(paste0("4.4.3_", prefix, "_cluster_PCA_tSNE_UMAP.pdf"), width = 60, height = 25)
plot_grid(ncol = 3,
  DimPlot(cluster_pscrna_seurat, reduction = "pca", group.by = "orig.ident"),
  DimPlot(cluster_pscrna_seurat, reduction = "tsne", group.by = "orig.ident", label = T),
  DimPlot(cluster_pscrna_seurat, reduction = "umap", group.by = "orig.ident", label = T),
  DimPlot(cluster_pscrna_seurat, reduction = "pca", group.by = "library_id"),
  DimPlot(cluster_pscrna_seurat, reduction = "tsne", group.by = "library_id", label = T),
  DimPlot(cluster_pscrna_seurat, reduction = "umap", group.by = "library_id", label = T),
  DimPlot(cluster_pscrna_seurat, reduction = "pca"),
  DimPlot(cluster_pscrna_seurat, reduction = "tsne", label = T),
  DimPlot(cluster_pscrna_seurat, reduction = "umap", label = T),
  DimPlot(cluster_pscrna_seurat, reduction = "pca", split.by = "orig.ident"),
  DimPlot(cluster_pscrna_seurat, reduction = "tsne", split.by = "orig.ident", label = T),
  DimPlot(cluster_pscrna_seurat, reduction = "umap", split.by = "orig.ident", label = T),
  DimPlot(cluster_pscrna_seurat, reduction = "pca", split.by = "library_id"),
  DimPlot(cluster_pscrna_seurat, reduction = "tsne", split.by = "library_id", label = T),
  DimPlot(cluster_pscrna_seurat, reduction = "umap", split.by = "library_id", label = T)
)
dev.off()
#=========================================================================================
# 4.5 Save
#=========================================================================================
save(cluster_pscrna_seurat, file=paste0(prefix, "_clusters_seurat.RData"), compress = F)
###########################################################################################
# 5. Identification of cluster markers gene
###########################################################################################
# Firstly, the expression of marker genes in the cluster was observed and visualized with vlnplot
## Visualization by cell clusters
for (i in 3:length(unique(mgene$name))) {
  cat(paste0("############### Plot VlnPlot for cluster ", unique(mgene$name)[i], " marker genes ###############\n"))
  try(VlnPlot(object = cluster_pscrna_seurat, features = filter(mgene, name == unique(mgene$name)[i])$gene, pt.size = 0, ncol = 5))
  ggsave(paste0("5_", prefix, "_cluster_", unique(mgene$name)[i], "_paper_markers_VlnPlot.pdf"), width = 30, height = 16, limitsize = FALSE)
}

# Visualization of all genes
paper_markers <- mgene %>% filter (! duplicated(genes)) %>% select(genes) %>% pull()

VlnPlot(object = cluster_pscrna_seurat, features = paper_markers, pt.size = 0, ncol = 5)
ggsave(paste0("5_", prefix, "_paper_markers_VlnPlot.pdf"), width = 50, height = 100, limitsize = FALSE)
#=========================================================================================
# 5.1 Identify all markers of each cluster
#=========================================================================================
#-----------------------------------------------------------------------------------------
# 5.1.1. Identify all markers of each cluster
#-----------------------------------------------------------------------------------------
# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(object = cluster_pscrna_seurat, 
                          only.pos = posType,
                          min.pct = pctmin,
                          thresh.use = usethresh,
                          logfc.threshold = thresholdlogfc)

dim(markers)
head(markers)
# The frequencies of all identified Merker genes were counted
table(table(markers$gene))
markers_annot <- markers %>% left_join(y = mgene, by = c("gene" = "genes")) %>% filter(!is.na(name))
markers_annot2 <- markers %>% left_join(y = unique(annotations[, c("gene_name", "Annotation_uniprot")]),
                                        by = c("gene" = "gene_name"))
write.csv(markers, paste0(prefix, "_Seurat_FindAllMarkers_each_cluster.csv"), row.names = FALSE)
write.csv(markers_annot, paste0(prefix, "_Seurat_FindAllMarkers_each_cluster_annot_markers.csv"), row.names = FALSE)
write.csv(markers_annot2, paste0(prefix, "_Seurat_FindAllMarkers_each_cluster_annot_uniprot.csv"), row.names = FALSE)
#-----------------------------------------------------------------------------------------
# 5.1.2. Get the marker gene that only appears once in all clusters
#-----------------------------------------------------------------------------------------
markers_single <- markers[markers$gene %in% names(table(markers$gene))[table(markers$gene) == 1],]
dim(markers_single)
table(table(markers_single$gene))
# The distribution of marker genes in each cluster was counted
table(markers_single$cluster)
head(markers_single)

topmaker <- function(n){
  topm <- markers_single %>% group_by(cluster) %>% top_n(n, avg_logFC)
  # create a scale.data slot for the selected genes
  topmSD <- ScaleData(cluster_pscrna_seurat, features = as.character(unique(topm$gene)), assay = "RNA")
  p1 <- DoHeatmap(object = cluster_pscrna_seurat, features = as.character(unique(topm$gene)))
  p11 <- DoHeatmap(object = topmSD, features = as.character(unique(topm$gene)))
  p2 <- VlnPlot(object = cluster_pscrna_seurat, features = as.character(unique(topm$gene)), pt.size = 0, ncol = 5)
  p22 <- VlnPlot(object = topmSD, features = as.character(unique(topm$gene)), pt.size = 0, ncol = 5)
  p3 <- FeaturePlot(cluster_pscrna_seurat, features = as.character(unique(topm$gene)), cols = c("lightgrey", "blue"), label = TRUE, ncol = 5)
  p33 <- FeaturePlot(topmSD, features = as.character(unique(topm$gene)), cols = c("lightgrey", "blue"), label = TRUE, ncol = 5)
  p4 <- DotPlot(cluster_pscrna_seurat, features = as.character(unique(topm$gene)), assay = "RNA") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  p44 <- DotPlot(topmSD, features = as.character(unique(topm$gene)), assay = "RNA") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  figures <- list(p1, p11, p2, p22, p3, p33, p4, p44)
  pdf(paste0("5.1.2_", prefix, "_Top", n, "_markers_DoHeatmap.pdf"), width = 30, height = 30)
  print(p1)
  print(p11)
  dev.off()
  pdf(paste0("5.1.2_", prefix, "_Top", n, "_markers_VlnPlot.pdf"), width = 50, height = 200)
  print(p2)
  print(p22)
  dev.off()
  pdf(paste0("5.1.2_", prefix, "_Top", n, "_markers_FeaturePlot.pdf"), width = 30, height = 30)
  print(p3)
  print(p33)
  dev.off()
  pdf(paste0("5.1.2_", prefix, "_Top", n, "_markers_DotPlot.pdf"), width = 30, height = 30)
  print(p4)
  print(p44)
  dev.off()
  pdf(paste0("5.1.2_", prefix, "_Top", n, "_markers_barplot_each_clusters.pdf"), width = 20, height = 10)
  mypar(3, 5, mar = c(4, 10, 3, 6))
  for (i in 1:length(unique(topm$cluster))) {
      barplot(sort(setNames(topm$avg_logFC, topm$gene)[topm$cluster == as.character(unique(topm$cluster))[i]], F), horiz = T, las = 1, main = paste0(as.character(unique(topm$cluster))[i], " vs. rest"), border = "white")
      abline(v = c(0, 0.25), lty = c(1, 2))
  }
  dev.off()
  return(figures)
}


for (n in c(topGenePlot)) {
  cat(paste0("########### Plot for Top", n, " genes ###########\n"))
  print(topmaker(n))
}
#-----------------------------------------------------------------------------------------
# 5.1.3. Annotation the average expression value of marker gene in the cluster and other clusters
#-----------------------------------------------------------------------------------------
# Get expression of genes for cells in and out of each cluster
getGeneClusterMeans <- function(gene, cluster){
  x <- GetAssayData(cluster_pscrna_seurat)[gene,]
  m <- tapply(x, ifelse(Idents(cluster_pscrna_seurat) == cluster, 1, 0), mean)
  mean.in.cluster <- m[2]
  mean.out.of.cluster <- m[1]
  return(list(mean.in.cluster = mean.in.cluster, mean.out.of.cluster = mean.out.of.cluster))
}

means <- mapply(getGeneClusterMeans, markers[,"gene"], markers[,"cluster"])
means <- matrix(unlist(means), ncol = 2, byrow = T)

colnames(means) <- c("mean.in.cluster", "mean.out.of.cluster")
rownames(means) <- markers[,"gene"]
markers2 <- cbind(markers, means)
head(markers2)
dim(markers2)
markers2_annot <- markers2 %>% left_join(y = mgene, by = c("gene" = "genes"))
write.csv(markers2_annot, paste0(prefix, "_Seurat_FindAllMarkers_each_cluster_expression.csv"), row.names = FALSE)
#-----------------------------------------------------------------------------------------
# 5.1.4. Do this for all different resolutions
#-----------------------------------------------------------------------------------------
markers_list <- list()
markers_annot <- list()
markers_annot2 <- list()
markers_single <- list()
topmaker <- list()
for (resolution in 1:length(resolutions)) {
  cat(paste0("########### Identify for resolution ", resolutions[resolution], " marker genes ###########\n"))
  Idents(cluster_pscrna_seurat) <- paste0("RNA_snn_res.", resolutions[resolution])
  paper_markers <- mgene %>% filter (! duplicated(genes), name != "BBM", name != "WUS") %>% select(genes) %>% pull()
  VlnPlot(object = cluster_pscrna_seurat, features = paper_markers, pt.size = 0, ncol = 5)
  ggsave(paste0("5.1.4_", prefix, "_resolution_", resolutions[resolution], "_paper_markers_VlnPlot.pdf"), width = 50, height = 100, limitsize = FALSE)
  markers_list[[resolution]] <- FindAllMarkers(object = cluster_pscrna_seurat, 
                                               only.pos = posType,
                                               min.pct = pctmin,
                                               thresh.use = usethresh,
                                               logfc.threshold = thresholdlogfc)

  markers_annot[[resolution]] <- markers_list[[resolution]] %>% left_join(y = mgene, by = c("gene" = "genes")) %>% filter(!is.na(name))
  markers_annot2[[resolution]] <- markers_list[[resolution]] %>% left_join(y = unique(annotations[, c("gene_name", "Annotation_uniprot")]), by = c("gene" = "gene_name"))
  write.csv(markers_list[[resolution]], paste0(prefix, "_resolution_", resolutions[resolution], "_Seurat_FindAllMarkers_each_cluster.csv"), row.names = FALSE)
  write.csv(markers_annot[[resolution]], paste0(prefix, "_resolution_", resolutions[resolution], "_Seurat_FindAllMarkers_each_cluster_annot_markers.csv"), row.names = FALSE)
  write.csv(markers_annot2[[resolution]], paste0(prefix, "_resolution_", resolutions[resolution], "_Seurat_FindAllMarkers_each_cluster_annot_uniprot.csv"), row.names = FALSE)

  markers_single[[resolution]] <- markers_list[[resolution]][markers_list[[resolution]]$gene %in% names(table(markers_list[[resolution]]$gene))[table(markers_list[[resolution]]$gene) == 1],]

  print(table(markers_single[[resolution]]$cluster))

  topmaker[[resolution]] <- function(n){

    topm <- markers_single[[resolution]] %>% group_by(cluster) %>% top_n(n, avg_logFC)

    # create a scale.data slot for the selected genes
    topmSD <- ScaleData(cluster_pscrna_seurat, features = as.character(unique(topm$gene)), assay = "RNA")
    p1 <- DoHeatmap(object = cluster_pscrna_seurat, features = as.character(unique(topm$gene)))
    p11 <- DoHeatmap(object = topmSD, features = as.character(unique(topm$gene)))
    p2 <- VlnPlot(object = cluster_pscrna_seurat, features = as.character(unique(topm$gene)), pt.size = 0, ncol = 5)
    p22 <- VlnPlot(object = topmSD, features = as.character(unique(topm$gene)), pt.size = 0, ncol = 5)
    p3 <- FeaturePlot(cluster_pscrna_seurat, features = as.character(unique(topm$gene)), cols = c("lightgrey", "blue"), label = TRUE, ncol = 5)
    p33 <- FeaturePlot(topmSD, features = as.character(unique(topm$gene)), cols = c("lightgrey", "blue"), label = TRUE, ncol = 5)
    p4 <- DotPlot(cluster_pscrna_seurat, features = as.character(unique(topm$gene)), assay = "RNA") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
    p44 <- DotPlot(topmSD, features = as.character(unique(topm$gene)), assay = "RNA") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
    figures <- list(p1, p11, p2, p22, p3, p33, p4, p44)
    cat(paste0("########### DoHeatmap ###########\n"))
    pdf(paste0("5.1.4_", prefix, "_resolution_", resolutions[resolution], "_Top", n, "_markers_DoHeatmap.pdf"), width = 30, height = 30)
    print(p1)
    print(p11)
    dev.off()
    cat(paste0("########### VlnPlot ###########\n"))
    pdf(paste0("5.1.4_", prefix, "_resolution_", resolutions[resolution], "_Top", n, "_markers_VlnPlot.pdf"), width = 50, height = 200)
    print(p2)
    print(p22)
    dev.off()
    cat(paste0("########### FeaturePlot ###########\n"))
    pdf(paste0("5.1.4_", prefix, "_resolution_", resolutions[resolution], "_Top", n, "_markers_FeaturePlot.pdf"), width = 30, height = 30)
    print(p3)
    print(p33)
    dev.off()
    cat(paste0("########### DotPlot ###########\n"))
    pdf(paste0("5.1.4_", prefix, "_resolution_", resolutions[resolution], "_Top", n, "_markers_DotPlot.pdf"), width = 30, height = 30)
    print(p4)
    print(p44)
    dev.off()
    return(figures)
  }


  for (n in c(topGenePlot)) {
    cat(paste0("########### Plot for Top", n, " and resolution ", resolutions[resolution], " genes ###########\n"))
    print(topmaker[[resolution]](n))
  }
}
Idents(cluster_pscrna_seurat) <- "RNA_snn_res.0.4"
#=========================================================================================
# 5.2 Identification of marker genes between specific clusters
#=========================================================================================
#-----------------------------------------------------------------------------------------
# 5.2.1. Default Wilcox method
#-----------------------------------------------------------------------------------------
# Determine differentiating markers for CD4+ T cell
selclustersw <- FindMarkers(cluster_pscrna_seurat,
                          ident.1 = 2,
                          ident.2 = c(0,4,10))
head(selclustersw)
# Add gene symbols to the DE table
selclustersw <- selclustersw %>%
              rownames_to_column(var = "gene") %>%
              left_join(y = unique(annotations[, c("gene_name", "Annotation_uniprot")]),
                        by = c("gene" = "gene_name"))

# Reorder columns and sort by padj      
selclustersw <- selclustersw[, c(1, 3:5,2,6:7)]
selclustersw <- selclustersw %>% dplyr::arrange(p_val_adj) 

# Candidate marker genes were selected for visualization
topFC <- selclustersw %>% top_n(topGenePlot, avg_logFC)
pdf(paste0("5.3.1_", prefix, "_Top",topGenePlot,"_wilcox_logFC_FindMarkers.pdf"), width = 20, height = 10)
for (gene in topFC$gene) {
  print(VlnPlot(cluster_pscrna_seurat, features = gene, pt.size = 0))
  print(FeaturePlot(cluster_pscrna_seurat, reduction = "umap", features = gene, label = TRUE, order = TRUE, min.cutoff = 'q10', repel = TRUE))
}
dev.off()

toppv <- selclustersw %>% top_n(topGenePlot, p_val_adj)
pdf(paste0("5.3.1_", prefix, "_Top",topGenePlot,"_wilcox_p_val_adj_FindMarkers.pdf"), width = 20, height = 10)
for (gene in toppv$gene) {
  print(VlnPlot(cluster_pscrna_seurat, features = gene, pt.size = 0))
  print(FeaturePlot(cluster_pscrna_seurat, reduction = "umap", features = gene, label = TRUE, order = TRUE, min.cutoff = 'q10', repel = TRUE))
}
dev.off()
#-----------------------------------------------------------------------------------------
# 5.2.2. roc method
#-----------------------------------------------------------------------------------------
selclusters <- FindMarkers(cluster_pscrna_seurat,
                           ident.1 = 2,
                           ident.2 = c(0,4,10),
                           test.use = "roc")

head(selclusters)
topr <- selclusters %>% top_n(topGenePlot, myAUC)

pdf(paste0("5.3.2_", prefix, "_Top",topGenePlot,"_Roc_FindMarkers.pdf"), width = 20, height = 10)
for (gene in row.names(topr)) {
  print(VlnPlot(cluster_pscrna_seurat, features = gene, pt.size = 0))
  print(FeaturePlot(cluster_pscrna_seurat, reduction = "umap", features = gene, label = TRUE, order = TRUE, min.cutoff = 'q10', repel = TRUE))
}
dev.off()
#=========================================================================================
# 5.3 Cluster modification
#=========================================================================================
#-----------------------------------------------------------------------------------------
# 5.3.1 Merge some clusters: based on TSNE and Heirarchical tree
#-----------------------------------------------------------------------------------------
# merge clusters 6 and 7 into 0 and cluster 9 into 13
experiment.merged <- RenameIdents(object = experiment.merged, '6' = '0', '7' = '0', '9' = '13')
# Look at the cell statistics for each cluster
table(Idents(experiment.merged))
DimPlot(object = experiment.examples, pt.size=0.5, label = T, reduction = "tsne")
#-----------------------------------------------------------------------------------------
# 5.3.2 Replace a cluster in this resolution with a cluster in another resolution (you can add some clusters in another resolution)
#-----------------------------------------------------------------------------------------
cluster_pscrna_seurat_man <- cluster_pscrna_seurat
newIdent <- as.character(Idents(cluster_pscrna_seurat_man))
newIdent[newIdent == '9'] <- paste0("R", as.character(cluster_pscrna_seurat_man$RNA_snn_res.0.8[newIdent == '9']))
Idents(cluster_pscrna_seurat_man) <- as.factor(newIdent)
table(Idents(cluster_pscrna_seurat_man))
#-----------------------------------------------------------------------------------------
# 5.3.3 Remove some clusters
#-----------------------------------------------------------------------------------------
# create a new tmp object with those removed
cluster_pscrna_seurat_man <- cluster_pscrna_seurat_man[,-which(Idents(cluster_pscrna_seurat_man) %in% 
                             c("14", "15", "R1", "R10", "R7", "R9"))]
table(Idents(cluster_pscrna_seurat_man))
table(Idents(cluster_pscrna_seurat))
#-----------------------------------------------------------------------------------------
# 5.3.4 Order clusters
#-----------------------------------------------------------------------------------------
Idents(cluster_pscrna_seurat_man) <- factor(cluster_pscrna_seurat_man@active.ident, 
                                            levels=c("0","1","2","3","4","5","6","7","8","R14","R18","10","11","12","13"))

levels(cluster_pscrna_seurat_man@active.ident)
#-----------------------------------------------------------------------------------------
# 5.3.5 tSNE and UMAP visualization of final cluster
#-----------------------------------------------------------------------------------------
umapcols <- length(table(Idents(cluster_pscrna_seurat_man)))
umapcolors <- colorRampPalette(brewer.pal(12, "Paired"), alpha=TRUE)(umapcols)
pdf(paste0("5.3.5_", prefix, "_Manual_order_cluster_UMAP.pdf"), width = 12, height = 10)
DimPlot(cluster_pscrna_seurat_man, reduction = "umap", pt.size = 1, label = TRUE, repel = T, label.box = T, label.size = 7, cols = umapcolors)
DimPlot(cluster_pscrna_seurat_man, pt.size=1, reduction = "umap", label = TRUE, repel = T, label.box = T, label.size = 7, cols = umapcolors, split.by = "orig.ident")
DimPlot(cluster_pscrna_seurat_man, pt.size=1, reduction = "umap", label = TRUE, repel = T, label.box = T, label.size = 7, cols = umapcolors, split.by = "library_id")
dev.off()

pdf(paste0("5.3.5_", prefix, "_Manual_order_cluster_UMAP2.pdf"), width = 20, height = 6)
DimPlot(cluster_pscrna_seurat_man, pt.size=1, reduction = "umap", label = F, cols = umapcolors, split.by = "orig.ident")
dev.off()

pdf(paste0("5.3.5_", prefix, "_Manual_order_cluster_UMAP3.pdf"), width = 40, height = 6)
DimPlot(cluster_pscrna_seurat_man, pt.size=1, reduction = "umap", label = F, cols = umapcolors, split.by = "library_id")
dev.off()

pdf(paste0("5.3.5_", prefix, "_Manual_order_cluster_tSNE.pdf"), width = 12, height = 10)
DimPlot(cluster_pscrna_seurat_man, reduction = "tsne", pt.size = 1, label = TRUE, repel = T, label.box = T, label.size = 7, cols = umapcolors)
DimPlot(cluster_pscrna_seurat_man, pt.size=1, reduction = "tsne", label = TRUE, repel = T, label.box = T, label.size = 7, cols = umapcolors, split.by = "orig.ident")
DimPlot(cluster_pscrna_seurat_man, pt.size=1, reduction = "tsne", label = F, cols = umapcolors, split.by = "orig.ident")
DimPlot(cluster_pscrna_seurat_man, pt.size=1, reduction = "tsne", label = TRUE, repel = T, label.box = T, label.size = 7, cols = umapcolors, split.by = "library_id")
DimPlot(cluster_pscrna_seurat_man, pt.size=1, reduction = "tsne", label = F, cols = umapcolors, split.by = "library_id")
dev.off()

pdf(paste0("5.3.5_", prefix, "_Manual_order_cluster_tSNE2.pdf"), width = 20, height = 6)
DimPlot(cluster_pscrna_seurat_man, pt.size=1, reduction = "tsne", label = F, cols = umapcolors, split.by = "orig.ident")
dev.off()

pdf(paste0("5.3.5_", prefix, "_Manual_order_cluster_tSNE3.pdf"), width = 40, height = 6)
DimPlot(cluster_pscrna_seurat_man, pt.size=1, reduction = "tsne", label = F, cols = umapcolors, split.by = "library_id")
dev.off()
#-----------------------------------------------------------------------------------------
# 5.3.6 paper marker visualization
#-----------------------------------------------------------------------------------------
VlnPlot(object = cluster_pscrna_seurat_man, features = paper_markers, pt.size = 0, ncol = 5)
ggsave(paste0("5.3.6_", prefix, "_Manual_paper_markers_VlnPlot.pdf"), width = 50, height = 100, limitsize = FALSE)
#-----------------------------------------------------------------------------------------
# 5.3.7 Save metadata after cluster operation
#-----------------------------------------------------------------------------------------
cluster_pscrna_seurat_man <- RenameIdents(object = cluster_pscrna_seurat_man, 
                                  "0" = "TM0",
                                  "1" = "TM1",
                                  "2" = "TM2",
                                  "3" = "TM3",
                                  "4" = "TM4",
                                  "5" = "TM5",
                                  "6" = "TM6",
                                  "7" = "TM7",
                                  "8" = "TM8",
                                  "R14" = "TM9",
                                  "10" = "TM10",
                                  "11" = "TM11",
                                  "12" = "TM12",
                                  "13" = "TM13",
                                  "R18" = "TM14")
save(cluster_pscrna_seurat_man, file=paste0(prefix, "_manual_clusters_seurat.RData"), compress = F)
#=========================================================================================
# 5.4 Rename all identities
#=========================================================================================
cluster_pscrna_seurat <- RenameIdents(object = cluster_pscrna_seurat, 
                                  "TM0" = "cluster_name0",
                                  "TM1" = "cluster_name1",
                                  "TM2" = "cluster_name2",
                                  "TM3" = "cluster_name3",
                                  "TM4" = "cluster_name4",
                                  "TM5" = "cluster_name5",
                                  "TM6" = "cluster_name6",
                                  "TM7" = "cluster_name7",
                                  "TM8" = "cluster_name8",
                                  "TM9" = "cluster_name9",
                                  "TM10" = "cluster_name10",
                                  "TM11" = "cluster_name11",
                                  "TM12" = "cluster_name12",
                                  "TM13" = "cluster_name13",
                                  "TM14" = "cluster_name14")
# Plot the UMAP
pdf(paste0("5.4_", prefix, "_rename_identities.pdf"), width = 15, height = 10)
DimPlot(object = cluster_pscrna_seurat, 
        reduction = "tsne", 
        label = TRUE,
        pt.size = 1,
        label.size = 3,
        repel = TRUE)

DimPlot(object = cluster_pscrna_seurat, 
        reduction = "umap", 
        label = TRUE,
        pt.size = 1,
        label.size = 3,
        repel = TRUE)

# Subsetting samples: observe the cluster distribution in a sample
cluster_pscrna_seurat_sample <- list()
for (i in unique(aggregation$orig.ident)) {
  cluster_pscrna_seurat_sample[[i]] <- subset(cluster_pscrna_seurat, orig.ident == i)
  p_tsne <- DimPlot(object = cluster_pscrna_seurat_sample[[i]], group.by = "RNA_snn_res.0.8", pt.size=0.5, label = TRUE, reduction = "tsne") + theme_classic() + ggtitle(paste0("Sample of ",i))
  print(p_tsne)
  p_umap <- DimPlot(object = cluster_pscrna_seurat_sample[[i]], group.by = "RNA_snn_res.0.8", pt.size=0.5, label = TRUE, reduction = "umap") + theme_classic() + ggtitle(paste0("Sample of ",i))
  print(p_umap)
}
dev.off()
#=========================================================================================
# 5.5 Save
#=========================================================================================
save(list=ls(), file=paste0(prefix, "_rename_markers_seurat.RData"), compress = F)
write_rds(cluster_pscrna_seurat, file = paste0(prefix, "_rename_markers_seurat.rds"), compress = "none")  
###########################################################################################
# Session information
###########################################################################################
sessionInfo()
