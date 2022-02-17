# Author:                      Job van Riet
# Date of  creation:           07-02-2020
# Function:                    Generating main figures of SRSR manuscript.

# Packages and general functions ------------------------------------------

library(Rsubread)
library(plyr)
library(dplyr)
library(DESeq2)
library(ggplot2)


# Create a GTF of the SRSR elements used in counting. ---------------------

# Convert into FeatureCounts annotation file.
SRSR.Anno <- data.frame(
  GeneID = data.SRSR.GRanges$CRISPR_Id,
  Chr = seqnames(data.SRSR.GRanges),
  Start = start(data.SRSR.GRanges),
  End = end(data.SRSR.GRanges),
  Strand = '.',
  stringsAsFactors=FALSE
)

SRSR.Anno <- SRSR.Anno[SRSR.Anno$Start != '.' & SRSR.Anno$End != '.',]

# Metadata of the NGS-ProToCoL samples.
sampleInfo <- readxl::read_xlsx('~/git/bitbucket/immunophenoprostate/misc/supTable1_SampleOverview.xlsx', sheet = 5, trim_ws = T)

# Input BAM files, keep only the correct samples.
files.bam <- list.files('/mnt/data/ccbc_environment/project/immunoProfiles_PCa/datasetNGSProtocol/bam_hg38/', full.names = T, pattern = 'sortedByCoord_markDup.bam$')
files.bam <- files.bam[grepl(paste(sampleInfo$Accession, collapse = '|'), files.bam)]
counts.SRSR <- Rsubread::featureCounts(files.bam, annot.ext=SRSR.Anno, primaryOnly = T, strandSpecific = 2, ignoreDup = T, isPairedEnd = T, nthreads = 40)
save(counts.SRSR, file = '/mnt/data/ccbc_environment/project/hCRISPR/externalData/NGSProtocol/counts.SRSR.RData')

# load('/mnt/data/ccbc_environment/project/hCRISPR/externalData/NGSProtocol/counts.SRSR.RData')

# Convert to matrix.
countMatrix <- counts.SRSR$counts
colnames(countMatrix) <- gsub('.Aligned.sortedByCoord.markDup.bam', '', colnames(countMatrix))
colnames(countMatrix) <- gsub('\\.', '-', colnames(countMatrix))

# Only keep included samples.
countMatrix <- countMatrix[,colnames(countMatrix) %in% sampleInfo$Accession]

# Perform DESeq2 on Tumor vs Normal.
DESeq2.SRSR.TvsN <- DESeq2::DESeqDataSetFromMatrix(countData = countMatrix, colData = sampleInfo[match(colnames(countMatrix), sampleInfo$Accession),], design = ~Type)
DESeq2.SRSR.TvsN <- DESeq2::DESeq(DESeq2.SRSR.TvsN, test = 'Wald', parallel = T)

# Get diff. results.
DESeq2.SRSR.results <- DESeq2::results(DESeq2.SRSR.TvsN, pAdjustMethod = 'BH', filterFun = IHW::ihw, parallel = T, tidy = F, contrast = c('Type', 'Tumor', 'NAP'))

# Shrink the LFC.
DESeq2.SRSR.results.LFC <- DESeq2::lfcShrink(DESeq2.SRSR.TvsN, res = DESeq2.SRSR.results, parallel = T,  coef = 'Type_Tumor_vs_NAP', type='normal')

# Add ENSEMBL as column.
DESeq2.SRSR.results.LFC$CRISPR_Id <- rownames(DESeq2.SRSR.results.LFC)

# Re-add the t-statistic
DESeq2.SRSR.results.LFC$stat <- DESeq2.SRSR.results[match(rownames(DESeq2.SRSR.results), rownames(DESeq2.SRSR.results.LFC)),]$stat

# Re-add the q weight.
DESeq2.SRSR.results.LFC$weight <- DESeq2.SRSR.results[match(rownames(DESeq2.SRSR.results), rownames(DESeq2.SRSR.results.LFC)),]$weight

# Add type of comparison
DESeq2.SRSR.results.LFC$coef <- 'Type_Tumor_vs_NAP'

# Convert to tibble.
DESeq2.SRSR.results.LFC <- tibble::as_tibble(DESeq2.SRSR.results.LFC)

# Check if statistically significant top candidate.
DESeq2.SRSR.results.LFC <- DESeq2.SRSR.results.LFC %>% dplyr::mutate(topDiffCandidate = ifelse(padj <= 0.05 & abs(log2FoldChange) >= .5 & baseMean >= 10, 'Yes', 'No'))

save(DESeq2.SRSR.results.LFC, file = '/mnt/data/ccbc_environment/project/hCRISPR/externalData/NGSProtocol/DESeq2.SRSR.results.LFC.RData')

# load('/mnt/data/ccbc_environment/project/hCRISPR/externalData/NGSProtocol/DESeq2.SRSR.results.LFC.RData')


# Volcano -----------------------------------------------------------------

data <- DESeq2.SRSR.results.LFC %>% 
  # Check direction.
  dplyr::mutate(Direction = 'Not signficant', 
                Direction = ifelse(log2FoldChange < 0, 'Down-regulated', 'Up-regulated'), 
                Direction = ifelse(topDiffCandidate == 'Yes', Direction, 'Not signficant'),
                log2FoldChange = as.numeric(log2FoldChange),
                padj = as.numeric(padj))

customColors <- ifelse(data$Direction == 'Not signficant', 'grey50', ifelse(data$Direction == 'Up-regulated', '#B02428', '#6697EA'))
names(customColors) <- data$Direction

EnhancedVolcano::EnhancedVolcano(data.frame(data), 
                                 lab = data$CRISPR_Id, 
                                 axisLabSize = 10, 
                                 title = NULL, subtitle = NULL, 
                                 FCcutoff = .5, pCutoff = 0.01,
                                 x = 'log2FoldChange', transcriptLabSize = 2, transcriptLabvjust = -.9, transcriptLabhjust = .5,
                                 y = 'padj', colCustom = customColors, colAlpha = 1,
                                 xlim = c(-3, 5),
                                 legendPosition = 'bottom', selectLab = data[data$topDiffCandidate == 'Yes',]$CRISPR_Id,
                                 legendLabSize = 10,
                                 legendIconSize = 3.0) +
  theme(      
    panel.grid.minor = element_line(colour = 'grey50', size = .3, linetype = 'dotted'),
    panel.grid.major = element_line(colour = 'grey70', size = .3, linetype = 'dotted')
  )

# Heatmap -----------------------------------------------------------------

DESeq2.SRSR.TvsN.vst <- DESeq2::vst(DESeq2.SRSR.TvsN[rowSums(DESeq2::counts(DESeq2.SRSR.TvsN, normalized = F)) > 1,], blind = T)

# Perform PCA on all genes.
normPCA <- prcomp(t(assay(DESeq2.SRSR.TvsN.vst)), center = T, scale. = T)

plot.all <- ggbiplot::ggbiplot(normPCA, choices = c(1,2), alpha = 0, groups = colData(DESeq2.SRSR.TvsN)$Type, var.axes = F, ellipse = T) +
  geom_point(aes(shape = colData(DESeq2.SRSR.TvsN)$Type, color = colData(DESeq2.SRSR.TvsN)$Type)) +
  scale_color_manual(name = NULL, values = c(NAP = 'darkblue', Tumor = '#ff6145')) +
  scale_shape_manual('Tissue', values = c(16,1)) +
  coord_cartesian() +
  theme(
    legend.position = 'bottom', 
    axis.title = element_text(size = 8),
    panel.grid.major.x = element_line(colour = 'grey90', linetype = 'longdash'),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'grey90', linetype = 'longdash'),
    panel.grid.minor.y = element_blank(),
    text=element_text(size=10, family='Helvetica'),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20')
  )


# Scree plot of first 20 PC
plotData <- factoextra::get_eigenvalue(normPCA); plotData$Dimension <- paste('PC', 1:nrow(plotData))
plot.Scree.All <- ggplot(plotData[1:20,], aes(x = reorder(Dimension, -variance.percent), y = variance.percent)) +
  geom_bar(stat = 'identity', position = 'dodge', width = .8, color= 'black', fill = 'tomato') +
  scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, 5), expand = c(.1,0)) +
  labs(x = 'Dimensions', y = 'Estimated percentage\nof explained variance') +
  geom_text(aes(label = round(..y.., 2)), stat='identity', position = position_dodge(width=0.8), vjust= .4, hjust = -.5,  angle = 90, size = 2.5) +
  theme(
    legend.position = 'none',
    text = element_text(size = 8, family='Helvetica', colour = 'black'),
    axis.text.x = element_text(angle = 90),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dotted'),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = NA, colour = NA)
  )


# PCA on only diff SRSR
z <- t(assay(DESeq2.SRSR.TvsN.vst[rownames(DESeq2.SRSR.TvsN.vst) %in% DESeq2.SRSR.results.LFC[DESeq2.SRSR.results.LFC$topDiffCandidate == 'Yes',]$CRISPR_Id,]))
colnames(z) <- sprintf('SRSR #%s', colnames(z))
normPCA.diff <- prcomp(z, center = T, scale. = T)

plot.diff <- ggbiplot::ggbiplot(normPCA.diff, choices = c(1,2), alpha = 0, groups = colData(DESeq2.SRSR.TvsN)$Type, var.axes = F, varname.size = 2, ellipse = T) +
  geom_point(aes(shape = colData(DESeq2.SRSR.TvsN)$Type, color = colData(DESeq2.SRSR.TvsN)$Type)) +
  scale_color_manual(name = NULL, values = c(NAP = 'darkblue', Tumor = '#ff6145')) +
  scale_shape_manual('Tissue', values = c(16,1)) +
  coord_cartesian() +
  theme(
    legend.position = 'bottom', 
    axis.title = element_text(size = 8),
    panel.grid.major.x = element_line(colour = 'grey90', linetype = 'longdash'),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'grey90', linetype = 'longdash'),
    panel.grid.minor.y = element_blank(),
    text=element_text(size=10, family='Helvetica'),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20')
  )


# Scree plot of first 20 PC
plotData <- factoextra::get_eigenvalue(normPCA.diff); plotData$Dimension <- paste('PC', 1:nrow(plotData))
plot.Scree.Diff <- ggplot(plotData[1:20,], aes(x = reorder(Dimension, -variance.percent), y = variance.percent)) +
  geom_bar(stat = 'identity', position = 'dodge', width = .8, color= 'black', fill = 'tomato') +
  scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, 5), expand = c(.1,0)) +
  labs(x = 'Dimensions', y = 'Estimated percentage\nof explained variance') +
  geom_text(aes(label = round(..y.., 2)), stat='identity', position = position_dodge(width=0.8), vjust= .4, hjust = -.5,  angle = 90, size = 2.5) +
  theme(
    legend.position = 'none',
    text = element_text(size = 8, family='Helvetica', colour = 'black'),
    axis.text.x = element_text(angle = 90),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dotted'),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = NA, colour = NA)
  )


cowplot::plot_grid(plot.all, plot.Scree.All, plot.diff, plot.Scree.Diff, align = 'hv', axis = 'tblr', nrow = 2, rel_widths = c(.45, .3))


# Heatmap -----------------------------------------------------------------

heatData <- assay(DESeq2::vst(DESeq2.SRSR.TvsN, blind = T))

# Show only the diff. SRSR.
heatData <- heatData[rownames(heatData) %in% DESeq2.SRSR.results.LFC[DESeq2.SRSR.results.LFC$topDiffCandidate == 'Yes',]$CRISPR_Id,]
rownames(heatData) <- sprintf('SRSR #%s', rownames(heatData))

# Heatmap colors.
paletteLength <- 101
heat.colors <- colorRampPalette(c('#6697EA', 'white', '#B02428'))(paletteLength)

# Column annotation.
annotation.col <- data.frame('Tissue' = factor(DESeq2.SRSR.TvsN$Type), row.names = DESeq2.SRSR.TvsN$Accession)
annotation.colors <- list('Tissue' = c(NAP = '#52639E', Tumor = '#EE6D7A'))

pheatmap::pheatmap(
  heatData, 
  fontsize_row = 4, treeheight_col = 15, treeheight_row = 15,
  angle_col = 45,  show_colnames = F, show_rownames = F,
  cluster_cols = T, annotation_col = annotation.col, annotation_colors = annotation.colors,
  cluster_rows = T, scale = 'row', cellheight = 2.5, cellwidth = 3, 
  clustering_method = 'ward.D2', 
  clustering_distance_cols = 'euclidean',
  clustering_distance_rows = 'euclidean', 
  border_color = 'grey90',
  color = heat.colors
)


# Plot some individual boxplots -------------------------------------------

boxData <- DESeq2::counts(DESeq2.SRSR.TvsN, normalized = T)

# Show only the diff. SRSR.
boxData <- boxData[rownames(boxData) %in% DESeq2.SRSR.results.LFC[DESeq2.SRSR.results.LFC$topDiffCandidate == 'Yes',]$CRISPR_Id,]

# Make a melted dataframe for ggplot2.
boxData <- reshape2::melt(boxData)
boxData <- boxData %>% dplyr::inner_join(tibble::as_tibble(SummarizedExperiment::colData(DESeq2.SRSR.TvsN)), by = c('Var2' = 'Accession'))

boxData %>% dplyr::filter(Var1 %in% head((DESeq2.SRSR.results.LFC %>% 
                                            dplyr::filter(topDiffCandidate == 'Yes') %>% 
                                            dplyr::filter(baseMean >= 20, lfcSE < 0.2, abs(log2FoldChange) >= 1) %>% 
                                            dplyr::arrange(padj) %>% dplyr::distinct(CRISPR_Id))$CRISPR_Id, n = 10)) %>% 
  dplyr::mutate(SRSR = sprintf('SRSR %s', Var1)) %>% 
  # Generate a boxplot.
  ggplot(., aes(SRSR, value, color = Type)) + 
  stat_boxplot(geom ='errorbar', position = position_dodge(.8), width = .25, alpha = .33) + 
  geom_boxplot(width = .5, position = position_dodge(.8), outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(aes(color = Type), dodge.width = .8, alpha = .2) +
  scale_color_manual( values = c(NAP = '#52639E', Tumor = '#EE6D7A'), guide = guide_legend(title = 'Tissue', title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 1, keyheight = 1)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(), limits = c(0, 500), breaks = c(0, 5, 10, 25, 50, 100, 250, 500, 750, 1000)) +
  labs(x = NULL, y = 'Normalized counts (log10)') + 
  theme(
    legend.position = 'bottom', 
    text = element_text(size = 8, family='Helvetica', colour = 'black'),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_blank(),
    axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
    panel.background = element_rect(fill = NA, colour = NA)
  )


# Write to table ----------------------------------------------------------

z <- DESeq2.SRSR.results.LFC %>% 
  dplyr::filter(topDiffCandidate == 'Yes') %>% 
  dplyr::group_by(CRISPR_Id) %>% 
  dplyr::summarise(info = sprintf('%s; %s', round(log2FoldChange, 2), gtools::stars.pval(padj))) %>% 
  dplyr::ungroup()

x <- readr::read_delim('srsr.txt', delim = '\t')
x %>% dplyr::left_join(z, by = c('SRSR ID' = 'CRISPR_Id')) %>% write.table(., file = 'asd.txt', quote = F, row.names = F, sep = '\t')
