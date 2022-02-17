# Author:                      Job van Riet
# Date of  creation:           05-02-2020
# Function:                    Import the MiTranscriptome into an RData object.


# Libraries ---------------------------------------------------------------

library(dplyr)
library(DESeq2)
library(Homo.sapiens)
library(ggplot2)

# Import ------------------------------------------------------------------

path.mitranscriptome.cnts <- '/mnt/data/ccbc_environment/project/hCRISPR/externalData/MiTranscriptome/mitranscriptome.cts.txt'
path.mitranscriptome.gtf <- '/mnt/data/ccbc_environment/project/hCRISPR/externalData/MiTranscriptome/mitranscriptome.v2.gtf'
path.mitranscriptome.samples <- '/mnt/data/ccbc_environment/project/hCRISPR/externalData/MiTranscriptome/samples.manifest.fix.tsv'


# Import Samples ----------------------------------------------------------

data.samples <- readr::read_delim(path.mitranscriptome.samples, delim = '\t', trim_ws = T, progress = T, col_names = T) %>% 
  dplyr::mutate(tissue = stringr::str_to_title(gsub('_', '-',tissue)), cancer_progression = stringr::str_to_title(gsub('_', ' ',cancer_progression))) %>%
  dplyr::filter(cancer_progression %in% c('Normal', 'Primary Cancer'))


# Clean MiTranscripts -----------------------------------------------------

# Read full GTF of all MiTranscripts and MiGenes.
data.gtf <- rtracklayer::import.gff(path.mitranscriptome.gtf)

# Make unique MiTranscript identifier.
data.gtf$MiTranscriptID <- paste(data.gtf$gene_id, data.gtf$tss_id, sep = '-')

# # Reduce transcript coordinates to gene-coordinates. (Very slow)
# cl <- BiocParallel::MulticoreParam(workers = 50, progressbar = T, stop.on.error = T)
# BiocParallel::register(cl, default = TRUE)
# 
# # Do the reduction per gene.
# idx <- base::split(seq_len(length(data.gtf)), list(data.gtf$gene_id))
# data.gtf.perGene <- GenomicRanges::GenomicRangesList(BiocParallel::bplapply(idx, function(x) IRanges::reduce(data.gtf[x,]), BPPARAM = cl))
# data.gtf.perGene.combined <- unlist(data.gtf.perGene)
# 
# # Only keep the gene-info of the genes in MiTranscriptome.
# data.gtf.perGene.combined <- data.gtf.perGene.combined[names(data.gtf.perGene.combined) %in% data.cnts$library_id,]
# data.gtf.perGene.combined <- data.gtf.perGene.combined[match(names(data.gtf.perGene.combined), data.cnts$library_id)]
# 
# mcols(data.gtf.perGene.combined)$MiID <- names(data.gtf.perGene.combined)
# 
# # Check if a MiTranscript overlaps with a known human gene.
# geneAnnotation <- rtracklayer::import.gff('/mnt/data/ccbc_environment/general/annotation/hg19/GENCODE/gencode.v28lift37.annotation.gtf.gz')
# geneAnnotation <- geneAnnotation[geneAnnotation$type == 'gene',]
# geneAnnotation <- geneAnnotation[geneAnnotation$level %in% c('1', '2'),]
# geneAnnotation <- geneAnnotation[!grepl('RP11-|CTD-', geneAnnotation$gene_name),]
# 
# overlapHumanGenes <- IRanges::findOverlaps(geneAnnotation, data.gtf.perGene.combined, type = 'any', minoverlap = 5)
# 
# z <- data.frame(SYMBOL = geneAnnotation[overlapHumanGenes@from]$gene_name, MiID = factor(names(data.gtf.perGene.combined[overlapHumanGenes@to])))
# z <- droplevels(z)
# z <- z %>% multidplyr::partition(MiID, cluster = multidplyr::create_cluster(30)) %>% dplyr::summarise(SYMBOL = paste(unique(SYMBOL), collapse = ', ')) %>% dplyr::collect()
# 
# mcols(data.gtf.perGene.combined) <- DataFrame(data.frame(mcols(data.gtf.perGene.combined)) %>% dplyr::left_join(z))
# 
# save(data.gtf.perGene.combined, file = '/mnt/data/ccbc_environment/project/hCRISPR/externalData/MiTranscriptome/data.gtf.perGene.combined.RData')

load('/mnt/data/ccbc_environment/project/hCRISPR/externalData/MiTranscriptome/data.gtf.perGene.combined.RData')

# Cohort selection --------------------------------------------------------

# Import raw counts of MiGenes (Huge file).
data.cnts <- readr::read_delim(path.mitranscriptome.cnts, delim = '\t', trim_ws = T, progress = T, col_names = T)

# Subset on samples for which we have metadata.
data.samples.subset <- data.samples[data.samples$library_id %in% colnames(data.cnts),]

# Only use cohorts for which we have enough tumor and normal samples. (>10)
sampleOverview <- data.samples.subset %>% dplyr::group_by(tissue) %>% 
  dplyr::summarise(Primary = sum(cancer_progression == 'Primary Cancer'),
                   Normal = sum(cancer_progression == 'Normal')) %>% 
  dplyr::filter(Normal >= 10) %>% reshape2::melt()

# Plot sample overview.
plot.sampleOverview <- ggplot(sampleOverview, aes(x = reorder(tissue, -value), y = value, fill = variable)) + geom_bar(stat = 'identity', color = 'grey25', width = .8, position = 'dodge') +
  scale_y_sqrt(name = '# Samples', breaks = c(0, 10, 25, 50, 100, 250, 500, 750, 1000, 1250), limit = c(0, 1250), expand = expand_scale(mult = c(0, .15))) +
  scale_fill_manual(values = c('tomato', 'skyblue'), labels = c('Malignant tissue', 'Healthy tissue')) +
  labs(x = 'MiTranscriptome cohorts') +
  geom_text(aes(label=value), position=position_dodge(width=0.8), vjust= -.2, hjust = -.2, angle = 45, size = 2.5) +
  guides( fill = guide_legend(title = 'Progression', title.position = 'top', title.hjust = 0.5, ncol = 2, keywidth = 0.5, keyheight = 0.5)) +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  theme(legend.position = 'bottom', axis.title.y = element_text(size = 11), text=element_text(size=10, family='Helvetica'),
        panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
        panel.grid.major.x = element_blank(),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        strip.background = element_rect(colour = 'grey20', fill = 'white')
  )


# Subset on SRSRs --------------------------------------------------------

# SRSR list (hg19) from script 1.

# Perform overlap.
data.gtf.perGene.combined.overlap <- IRanges::subsetByOverlaps(data.gtf.perGene.combined, data.SRSR.GRanges.hg19, minoverlap = 5, ignore.strand = T)
data.cnts <- data.cnts[data.cnts$library_id %in% names(data.gtf.perGene.combined.overlap),]

# Perform DESeq per tissue ------------------------------------------------

data.samples.subset.minNormals <- data.samples.subset %>% filter(tissue %in% (sampleOverview %>% dplyr::filter(variable == 'Normal', value >= 10))$tissue)
data.samples.subset.minNormals <- data.samples.subset.minNormals %>% dplyr::mutate(cancer_progression = gsub(' ', '_', cancer_progression))

# Create cluster to perform in parallel (if requested)
cl <- BiocParallel::MulticoreParam(workers = 3, progressbar = T)
BiocParallel::register(BiocParallel::bpstart(cl))

# Run per cohort.
allDESeq.MiTranscriptome <- BiocParallel::bplapply(unique(data.samples.subset.minNormals$tissue), function(tissue){

  # Get samples per tissue.
  data.samples.subset.minNormals.tissue <- data.samples.subset.minNormals[data.samples.subset.minNormals$tissue == tissue,]

  data.cnts.subset.tissue <- data.cnts[colnames(data.cnts) %in%  c('library_id', data.samples.subset.minNormals.tissue$library_id)]
  data.samples.subset.minNormals.tissue <- data.samples.subset.minNormals.tissue[match(colnames(data.cnts.subset.tissue), data.samples.subset.minNormals.tissue$library_id),]
  data.samples.subset.minNormals.tissue <- data.samples.subset.minNormals.tissue[complete.cases(data.samples.subset.minNormals.tissue),]

  # Subset samples.
  data.cnts.subset.tissue <- data.frame(data.cnts.subset.tissue, check.names = F)
  rownames(data.cnts.subset.tissue) <- data.cnts.subset.tissue$library_id
  data.cnts.subset.tissue$library_id <- NULL

  # Generate DESeq object.
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = data.cnts.subset.tissue, colData = data.samples.subset.minNormals.tissue, design = ~cancer_progression, tidy = F)

  # Analysis
  dds <- DESeq2::DESeq(dds, parallel = T, quiet = F, BPPARAM = BiocParallel::MulticoreParam(12))

  return(dds)
}, BPPARAM = cl)

# Terminate cluster.
BiocParallel::bpstop(cl)

# Correct names.
names(allDESeq.MiTranscriptome) <- unique(data.samples.subset.minNormals$tissue)

# Save image.
save(allDESeq.MiTranscriptome, file = '/mnt/data/ccbc_environment/project/hCRISPR/externalData/MiTranscriptome/DESeqPerMiTranscriptomeTissue.RData')
