# Author:                      Job van Riet
# Date of  creation:           05-02-2020
# Function:                    Generating MITranscriptome figures of SRSR manuscript.


# Packages and general functions ------------------------------------------

library(plyr)
library(dplyr)
library(ggplot2)
library(DESeq2)


# Import MiTranscriptome (DESeq) and SRSRs -----------------------------

# Load DESeq object.
load('/mnt/data/ccbc_environment/project/hCRISPR/externalData/MiTranscriptome/DESeqPerMiTranscriptomeTissue.RData')

# Load MiTranscriptome transcript information.
load('/mnt/data/ccbc_environment/project/hCRISPR/externalData/MiTranscriptome/data.gtf.perGene.combined.RData')


# Filter MiTranscriptome on diff. transcripts -----------------------------

results.MiTranscriptome <- lapply(names(allDESeq.MiTranscriptome), function(tissue){
  test <- allDESeq.MiTranscriptome[[tissue]]
  dds.results <- tibble::as_tibble(DESeq2::results(test, pAdjustMethod = 'BH', tidy = T, contrast = c('cancer_progression', 'Primary_Cancer', 'Normal')))
  dds.results <- dds.results %>% 
    dplyr::mutate(
      Direction = ifelse(log2FoldChange < 0, 'Down-regulated in cancer', 'Up-regulated in cancer'), 
      tissue = tissue,
      MiID = row,
      row = NULL)
  return(dds.results)
})

results.MiTranscriptome <- do.call(rbind, results.MiTranscriptome)

# Add the overlapping CRISPR identifier.
overlapMiTranscriptome <- tibble::as_tibble(IRanges::findOverlaps(data.gtf.perGene.combined, data.SRSR.GRanges.hg19, minoverlap = 5, ignore.strand = T))
overlapMiTranscriptome$MiID <- data.gtf.perGene.combined[overlapMiTranscriptome$queryHits,]$MiID
overlapMiTranscriptome$SYMBOL <- data.gtf.perGene.combined[overlapMiTranscriptome$queryHits,]$SYMBOL
overlapMiTranscriptome$CRISPR_Id <- data.SRSR.GRanges.hg19[overlapMiTranscriptome$subjectHits]$CRISPR_Id

overlapMiTranscriptome <- overlapMiTranscriptome %>% dplyr::group_by(MiID) %>% dplyr::summarise(
  CRISPR_Id = paste(unique(CRISPR_Id), collapse = ', '),
  Identifier = sprintf('%s (%s)', unique(MiID), ifelse(is.na(unique(SYMBOL)), 'Unannotated', unique(SYMBOL)))
)

# Combine data.
results.MiTranscriptome <- results.MiTranscriptome %>% dplyr::left_join(overlapMiTranscriptome, by = c('MiID' = 'MiID'))

# Filter on diff. expression.
results.MiTranscriptome.filtered <- results.MiTranscriptome %>% dplyr::filter(padj <= 0.01, baseMean >= 10, abs(log2FoldChange) >= 1)


# Overview of diff. MiTranscript per tissue -------------------------------

plotData <- results.MiTranscriptome.filtered %>% dplyr::group_by(tissue, Direction) %>% dplyr::tally() %>% dplyr::ungroup()
plotData <- plotData %>% dplyr::mutate(
  tissue = factor(tissue, levels = (plotData %>% dplyr::filter(grepl('Up-reg', Direction)) %>% dplyr::arrange(-n))$tissue),
  Direction = factor(Direction, levels = c('Up-regulated in cancer', 'Down-regulated in cancer'))
)

ggplot(plotData, aes(x = tissue, y = n, label = n, fill = tissue, group = Direction)) + 
  geom_bar(stat = 'identity', col = 'black', width = .75, position = position_dodge()) +
  scale_fill_manual(values = unname(hues::iwanthue(12)), guide = F) +
  scale_y_continuous( limits = c(0, 1500), breaks = c(0, 250, 500, 750, 1000, 1250, 1500)) + 
  geom_text(stat= 'identity', position = position_dodge(width = .75), vjust=-.75, size = 2.5) +
  labs(y = '# MiTranscripts overlapping with SRSR elements\n(q ≤ 0.01; Avg. Mean ≥ 10; logFC ≥ 1)') + 
  theme(
    legend.position = 'bottom',
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 7),
    text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'grey60', linetype = 'dashed'),
    panel.grid.minor.y = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )


# Heatmap of logFC per tissue per MiTranscript ----------------------------

heatData <- reshape2::dcast(results.MiTranscriptome %>% dplyr::select(log2FoldChange, tissue, Identifier), Identifier ~ tissue, value.var = 'log2FoldChange', fun.aggregate = sum)
heatData <- heatData[complete.cases(heatData),]
rownames(heatData) <- heatData$Identifier; heatData$Identifier <- NULL
heatData[is.na(heatData)] <- 0

# Filter on diff. MiTranscripts.
heatData <- heatData[rownames(heatData) %in% (results.MiTranscriptome.filtered)$Identifier,]

# Make upper/lower limit of logFC.
heatData[heatData > 4] <- 4.5
heatData[heatData < -4] <- -4.5

# Get the highest number of the matrix to create an even color scale.
paletteLength <- 100
heat.colors <- colorRampPalette(c('#4CAF50', 'lightgreen', 'white', 'tomato', '#F44336'))(101)

# Plot.
pheatmap::pheatmap(
  heatData, clustering_method = 'ward.D2', clustering_distance_cols = 'euclidean', clustering_distance_rows = 'euclidean',
  scale = 'none',
  treeheight_col = 25, treeheight_row = 25, 
  #cellheight = 7, cellwidth = 7, fontsize = 7,
  color = heat.colors, drop_levels = T, 
  border_color = 'grey75',
  show_colnames = T, 
  show_rownames = F,
  cluster_cols = T, 
  cluster_rows = T
)

# Tissue-specific markers -------------------------------------------------

dataBox <- results.MiTranscriptome.filtered %>% dplyr::arrange(abs(log2FoldChange)) %>% dplyr::group_by(tissue) %>% dplyr::filter(grepl('Unann', Identifier),  !grepl(',', CRISPR_Id), baseMean >= 75, abs(log2FoldChange) >= 2) %>% dplyr::filter(row_number() %in% 1:4) %>% dplyr::ungroup()

# Identifier == G008259 (Unannotated) & grepl(Lung, tissue)
# Head and Neck: G089948 (Unannotated) | G018828 (Unannotated)
# Thyroid: G046073 (Unannotated)
# Liver: G037026 (Unannotated)
# Uterus: G081524 (Unannotated) | G028957 (Unannotated)
# Stomach: G050669 (Unannotated)
# Colorectal: G013027 (Unannotated)
# Breast: G006354 (Unannotated)
# Bladder: G005582 (Unannotated)
# Kidney: G048595 (Unannotated)
# Esophagus: G005582 (Unannotated)
# Prostate: G052336 (Unannotated)
# Bladder: G009649 (Unannotated)

dataBox <- results.MiTranscriptome.filtered %>% dplyr::filter(
  (Identifier == 'G008259 (Unannotated)' & grepl('Lung', tissue)) |
    (Identifier == 'G083361 (Unannotated)' & grepl('Head', tissue)) |
    (Identifier == 'G046073 (Unannotated)' & grepl('Thyroid', tissue)) |
    (Identifier == 'G007433 (Unannotated)' & grepl('Liver', tissue)) |
    (Identifier == 'G007433 (Unannotated)' & grepl('Uterus', tissue)) |
    (Identifier == 'G050669 (Unannotated)' & grepl('Stomach', tissue)) |
    (Identifier == 'G013027 (Unannotated)' & grepl('Colorectal', tissue)) |
    (Identifier == 'G086960 (Unannotated)' & grepl('Breast', tissue)) |
    (Identifier == 'G048595 (Unannotated)' & grepl('Kidney', tissue)) |
    (Identifier == 'G005582 (Unannotated)' & grepl('Esophagus', tissue)) |
    (Identifier == 'G052336 (Unannotated)' & grepl('Prostate', tissue)) |
    (Identifier == 'G009649 (Unannotated)' & grepl('Bladder', tissue))
)

x <- apply(dataBox, 1, function(x){
  plotCounts(x['MiID'], x['tissue'], x['CRISPR_Id'], round(as.numeric(x['log2FoldChange']), 2), formatC(as.numeric(x['padj']), format = 'e', digits = 1), x['Identifier'])
})

plotCounts <- function(geneID, tissue, SRSR, logFC, q, identifier){
  normalizedCounts <- DESeq2::counts(allDESeq.MiTranscriptome[names(allDESeq.MiTranscriptome) == tissue][[1]], normalized = T)
  normalizedCounts <- normalizedCounts[rownames(normalizedCounts) == geneID,]
  normalizedCounts <- reshape2::melt(normalizedCounts)
  
  colData <- colData(allDESeq.MiTranscriptome[names(allDESeq.MiTranscriptome) == tissue][[1]])[c('library_id', 'cancer_progression')]
  normalizedCounts <- merge(normalizedCounts, colData, by.x = 0, by.y = 'library_id')
  normalizedCounts <- tibble::as_tibble(normalizedCounts)
  normalizedCounts <- normalizedCounts %>% dplyr::mutate(cancer_progression = factor(cancer_progression, levels = c('Primary_Cancer', 'Normal')))
  
  ggplot(normalizedCounts, aes(x = cancer_progression, y = value)) + 
    ggbeeswarm::geom_quasirandom(aes(color = cancer_progression), dodge.width = .8, alpha = .3) +
    stat_boxplot(geom ='errorbar', position = position_dodge(.8), width = .25, alpha = .33) + 
    geom_boxplot(width = .4, alpha = .33, position = position_dodge(.8), outlier.shape = NA) +
    scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 10, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000, 200000), limits = c(0, 400000)) +
    ggsignif::stat_signif(comparisons = list(c('Normal', 'Primary_Cancer')), test = 'wilcox.test', map_signif_level=TRUE) +
    scale_color_manual(labels = c(sprintf('Primary Cancer\n(%s)', tissue), sprintf('Healthy tissue\n(%s)', tissue)), values = c('tomato', 'skyblue'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.75, keyheight = 0.75)) +
    labs(x = sprintf('%s\nSRSR: #%s\nlogFC: %s; q: %s', identifier, SRSR, logFC, q), y = 'Normalized counts\n(RNA Expression)') + 
    theme(legend.position = 'bottom', 
          axis.title.y = element_text(size = 7), 
          axis.title.x = element_text(size = 7), 
          text=element_text(size = 7, family='Helvetica'),
          panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
          panel.grid.major.x = element_blank(),
          axis.text.x= element_blank(),
          axis.ticks.x=element_blank(),
          panel.background = element_rect(fill = 'white', colour = NA),
          panel.border = element_rect(fill = NA, colour = 'grey20')
    )
}

cowplot::plot_grid(plotlist = x, ncol = 6, align = 'hv', axis = 'tblr')


# Write table -------------------------------------------------------------

z <- results.MiTranscriptome %>% tidyr::separate_rows(CRISPR_Id, sep = ", ") %>% 
  dplyr::select(Identifier, CRISPR_Id, tissue, log2FoldChange, padj) %>% 
  dplyr::filter(padj <= 0.05, abs(log2FoldChange) >= 0.5) %>% 
  dplyr::group_by(CRISPR_Id, tissue) %>% 
  dplyr::summarise(info = sprintf('%s (%s; %s; %s)', unique(tissue), gsub(' \\(.*', '', Identifier), round(log2FoldChange, 2), gtools::stars.pval(padj))) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(CRISPR_Id) %>% 
  dplyr::summarise(info2 = paste(unique(info), collapse = ', ')) %>% 
  dplyr::ungroup()

x <- readr::read_delim('srsr.txt', delim = '\t')
x %>% dplyr::left_join(z, by = c('SRSR ID' = 'CRISPR_Id')) %>% write.table(., file = 'asd.txt', quote = F, row.names = F, sep = '\t')

