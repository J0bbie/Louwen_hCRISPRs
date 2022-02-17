# Author:                      Job van Riet
# Date of  creation:           03-02-2020
# Function:                    Generating main figures of SRSR manuscript.


# Packages and general functions ------------------------------------------

# Load libraries.
library(plyr)
library(dplyr)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Function to get overlap of SRSR with given GRanges.
findOverlapSRSR <- function(SRSR, data, minoverlap = 5){
  x <- GenomicRanges::findOverlaps(SRSR, data, minoverlap = minoverlap)
  y <- data[x@to]
  y$CRISPR_Id <- SRSR[x@from]$CRISPR_Id
  
  return(y)
}


# Import data -------------------------------------------------------------

# Import hg38 CRISPRCASFinder results.
data.SRSR <- readr::read_delim('/mnt/data/ccbc_environment/project/hCRISPR/CRISPRCasFinder/hg38_2019_11_29/TSV/Crisprs_REPORT.tsv', delim = '\t', trim_ws = T)

# Convert to GRanges.
data.SRSR.GRanges <- GenomicRanges::GRanges(seqnames = data.SRSR$Sequence, ranges = IRanges(data.SRSR$CRISPR_Start, data.SRSR$CRISPR_End))
GenomicRanges::mcols(data.SRSR.GRanges) <- data.SRSR
seqlevelsStyle(data.SRSR.GRanges) = "UCSC"

# Remove haplo-chromosomes, if any.
data.SRSR.GRanges <- R2CCBC::cleanSeqlevels(GenomicRanges::sort(data.SRSR.GRanges), genome ='hg38')


# Liftover the SRSR elements ----------------------------------------------

library(rtracklayer)

# hg38 to hg19.
ch <- rtracklayer::import.chain('/mnt/data/ccbc_environment/general/annotation/liftOversChains/hg38ToHg19.over.chain')
data.SRSR.GRanges.hg19 <- unlist(rtracklayer::liftOver(data.SRSR.GRanges, ch))
data.SRSR.GRanges.hg19$hg19Location <- sprintf('%s:%s-%s', seqnames(data.SRSR.GRanges.hg19), start(data.SRSR.GRanges.hg19), end(data.SRSR.GRanges.hg19))
data.SRSR.GRanges.hg19 <- data.SRSR.GRanges.hg19[!duplicated(data.SRSR.GRanges.hg19$CRISPR_Id)]

# Add the hg19 region back to the original list.
data.SRSR.GRanges$hg19Location <- NA
data.SRSR.GRanges[match(data.SRSR.GRanges.hg19$CRISPR_Id, data.SRSR.GRanges$CRISPR_Id),]$hg19Location <- data.SRSR.GRanges.hg19$hg19Location


# Overlap with known genes ------------------------------------------------

# Load gene information, only take genes with validated or manuallu-curated support.
geneAnnotation <- rtracklayer::import.gff('/mnt/data/ccbc_environment/general/annotation/hg38/GENCODE/gencode.v33.annotation.gtf')
geneAnnotation <- geneAnnotation[geneAnnotation$type == 'gene',]
geneAnnotation <- geneAnnotation[geneAnnotation$level %in% c('1', '2'),]
geneAnnotation <- geneAnnotation[geneAnnotation$gene_type != 'TEC',]


# Figure 2 - Location of SRSR elements ------------------------------------

# Bin hCRISPRs per 1Mb
SRSR.Binned <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(data.SRSR.GRanges), tilewidth=1E6, cut.last.tile.in.chrom=TRUE)
SRSR.Binned$totalCount <- IRanges::countOverlaps(SRSR.Binned, data.SRSR.GRanges)

SRSR.Binned[SRSR.Binned$totalCount > 15,]$totalCount <- 16

# Convert to data.frame in order to plot in Circos
locData <- data.frame(chr = GenomeInfoDb::seqnames(SRSR.Binned), start = IRanges::start(SRSR.Binned), end = IRanges::end(SRSR.Binned), count = SRSR.Binned$totalCount)
locData <- locData[locData$count > 0,]

# Circos ------------------------------------------------------------------

circlize::circos.clear()

# Start with a small gap after the last chromosome to allow for y-axis notations.
circlize::circos.par(start.degree = 90, gap.after = c(rep(1, length(GenomeInfoDb::seqlevels(data.SRSR.GRanges)) - 1), 10))
circlize::circos.initializeWithIdeogram(ideogram.height = 0.03, species = 'hg38', chromosome.index = GenomeInfoDb::seqlevels(data.SRSR.GRanges))

# Draw overlapping SRSR per bin.
circlize::circos.genomicTrackPlotRegion(locData, numeric.column = 4, ylim = c(0, 18), panel.fun = function(region, value, ...) {
  
  # Add segments as black points.
  circlize::circos.genomicRect(region, value[[1]], ytop = value[[1]], ybottom = 0, col = ifelse(value[[1]] > 15, '#F69F9A', '#101255'), border = ifelse(value[[1]] > 15, '#F69F9A', '#101255'))
  
  # Add axis lines.
  circlize::circos.lines(circlize::CELL_META$cell.xlim, c(0, 0), lty = 2, col = '#00000040')
  circlize::circos.lines(circlize::CELL_META$cell.xlim, c(5, 5), lty = 2, col = '#00000040')
  circlize::circos.lines(circlize::CELL_META$cell.xlim, c(10, 10), lty = 2, col = '#00000040')
  circlize::circos.lines(circlize::CELL_META$cell.xlim, c(15, 15), lty = 2, col = '#00000040')
  
}, track.height = circlize::uh(20))

circlize::circos.yaxis(side = 'left', at = c(0, 5, 10, 15), labels = c(0, 5, 10, '≥15'), sector.index = circlize::get.all.sector.index()[1], labels.cex = 0.3)


# Overlap of SRSR with genes ----------------------------------------------

# Overlap of hg38 genes.
overlapHumanGenes <- IRanges::subsetByOverlaps(geneAnnotation, data.SRSR.GRanges, type = 'any', minoverlap = 5)

# Non-overlap of SRSR with genes.
noHit <- length(data.SRSR.GRanges) - length(IRanges::subsetByOverlaps(data.SRSR.GRanges, geneAnnotation, invert = T, type = 'any', minoverlap = 5))

# Count frequency of gene-type.
overlapHumanGenesPerType <- tibble::as_tibble(data.frame(table(overlapHumanGenes$gene_type)))

# Summarize < 5 genes per category into the 'Other' group.
overlapHumanGenesPerType <- overlapHumanGenesPerType %>% dplyr::group_by(Var1) %>% dplyr::mutate(Type = ifelse(sum(Freq) < 5, 'Other', as.character(Var1)),
                                                                                                 Type = gsub('_', ' ', R.utils::capitalize(Type)),
                                                                                                 Type = ifelse(grepl('pseudo', Type), 'Pseudogene', Type)) %>% dplyr::ungroup()

overlapHumanGenesPerType <- overlapHumanGenesPerType %>% dplyr::group_by(Type) %>% dplyr::summarise(Freq = sum(as.numeric(Freq))) %>% dplyr::ungroup()

# Add SRSR without overlap.
overlapHumanGenesPerType <- rbind(overlapHumanGenesPerType, c('No overlap', noHit)) %>% 
  dplyr::mutate(Freq = as.numeric(Freq))

plot.geneCat <- ggplot(overlapHumanGenesPerType, aes(x = reorder(Type, Freq), y = Freq, label = Freq)) + 
  geom_bar(stat = 'identity', color = 'black', fill = '#005288', width = .8) +
  scale_y_continuous(limits = c(0, 15000), trans = scales::pseudo_log_trans(), breaks = c(0, 10, 50, 100, 250, 500, 1000, 2500, 5000, 10000)) + 
  labs(x = 'Categories', y = 'Frequency\n(log10)') +
  geom_text(color='black', size = 2.5, vjust = -1) +
  geom_abline(slope = 0) + 
  theme(
    legend.position = 'bottom', 
    axis.title.y = element_text(size = 8), 
    text = element_text(size=10, family='Helvetica'), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_line(linetype = 'dotted'),
    panel.background = element_rect(fill = NA, colour = NA), 
    panel.border = element_rect(fill = NA, colour = NA),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )


# Location of SRSR near or within genes -----------------------------------

location.withinGene <- ChIPseeker::annotatePeak(unique(subsetByOverlaps(data.SRSR.GRanges, geneAnnotation, type = 'any', minoverlap = 5)), tssRegion=c(-0, 0), level = 'gene', TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db')

overlapLocs <- tibble::as_tibble(data.frame(location.withinGene@annoStat)) %>% 
  dplyr::mutate(Frequency.Abs = round(Frequency / 100 * location.withinGene@peakNum)
  )

plot.geneLocs <- ggplot(overlapLocs, aes(x = reorder(Feature, Frequency.Abs), y = Frequency.Abs, label = Frequency.Abs)) + 
  geom_bar(stat = 'identity', color = 'black', fill = '#005288', width = .8) +
  scale_y_continuous(limits = c(0, 15000), trans = scales::pseudo_log_trans(), breaks = c(0, 10, 50, 100, 250, 500, 1000, 2500, 5000, 10000)) + 
  labs(x = 'Features', y = 'Frequency\n(log10)') +
  geom_text(color='black', size = 2.5, vjust = -1) +
  geom_abline(slope = 0) + 
  theme(
    legend.position = 'bottom', 
    axis.title.y = element_text(size = 8), 
    text = element_text(size=10, family='Helvetica'), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_line(linetype = 'dotted'),
    panel.background = element_rect(fill = NA, colour = NA), 
    panel.border = element_rect(fill = NA, colour = NA),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )


# Overlaps human repeats with the CRISPR-like elements --------------------

# Repeats from RepeatMasker (From NCBI: GCF_000001405.39)
humanRepeats <- readr::read_table('/mnt/data/ccbc_environment/project/hCRISPR/externalData/GCF_000001405.39_GRCh38.p13_rm.out', skip = 2, 
                                  col_names = c('SWScore', 'percDiv', 'percDel', 'percIns', 'chr', 'start', 'end', 'queryLeft', 'matchingRepeat', 'repeatFamily', 'startRepeat', 'endRepeat', 'ID'))

# Convert to GRanges and only keep the standard chromosomes.
humanRepeats <- GenomicRanges::makeGRangesFromDataFrame(humanRepeats, keep.extra.columns = T)
humanRepeats <- dropSeqlevels(humanRepeats, seqlevels(humanRepeats)[!grepl('NC_', seqlevels(humanRepeats))], pruning.mode = 'coarse')
seqlevels(humanRepeats) <- gsub('\\..*', '', gsub('NC_0*0', '', seqlevels(humanRepeats)))

# Convert NCBI seqlevels to chromosomes names.
seqlevelsStyle(humanRepeats) <- 'UCSC'
seqlevels(humanRepeats)[23:24] <- c('chrX', 'chrY')

# Summarize total repeats in human genome.
humanRepeats.total <- tibble::as_tibble(mcols(humanRepeats)) %>% dplyr::group_by(repeatFamily) %>% dplyr::summarise(repeatFamilyTotalFreq = n(), repeatFamilyTotal = sprintf('%s (%s)', unique(repeatFamily), n()))

# Overlap (min. 5 bp)
humanRepeats.overlap <- IRanges::findOverlaps(data.SRSR.GRanges, humanRepeats, minoverlap = 5)
humanRepeats.overlap <- tibble::as_tibble(humanRepeats.overlap)
humanRepeats.overlap$repeatFamily <- humanRepeats[humanRepeats.overlap$subjectHits]$repeatFamily

# Determine overlap / non-overlap
totalRepeatOverlap <- data.frame(
  Overlapping = length(IRanges::subsetByOverlaps(data.SRSR.GRanges, humanRepeats, minoverlap = 5)),
  'Non-overlapping' = length(IRanges::subsetByOverlaps(data.SRSR.GRanges, humanRepeats, invert = T, minoverlap = 5))
)

plot.repeats <- ggplot(reshape2::melt(totalRepeatOverlap), aes(x = reorder(variable, -value), y = value, label = value)) + 
  geom_bar(stat = 'identity', color = 'black', fill = '#005288', width = .8) +
  scale_y_continuous(limits = c(0, 10000)) + 
  labs(x = 'Overlap of SRSR\nwith RepeatMasker', y = 'Frequency') +
  geom_text(color='black', size = 2.5, vjust = -1) +
  geom_abline(slope = 0) + 
  theme(
    legend.position = 'bottom', 
    axis.title.y = element_text(size = 8), 
    text = element_text(size=10, family='Helvetica'),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = NA, colour = NA), 
    panel.border = element_rect(fill = NA, colour = NA),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )

# Convert data for plotting.
plotData <- humanRepeats.overlap %>% 
  dplyr::group_by(repeatFamily) %>% 
  dplyr::summarise(totalRepeat = n()) %>% 
  dplyr::ungroup()

# Add total repeat elements to name.
plotData <- plotData %>% dplyr::inner_join(humanRepeats.total)

plot.repeatsOverview <- ggplot(plotData %>% dplyr::filter(totalRepeat >= 75), aes(x = reorder(repeatFamily, -totalRepeat), y = totalRepeat, label = totalRepeat)) + 
  geom_bar(stat = 'identity', color = 'black', fill = '#005288', width = .8) +
  scale_y_continuous(limits = c(0, 600)) + 
  labs(x = 'Repeat Family', y = 'Frequency') +
  geom_text(color='black', size = 2.5, vjust = -1) +
  geom_abline(slope = 0) + 
  theme(
    legend.position = 'bottom', 
    axis.title.y = element_text(size = 8), 
    text = element_text(size=10, family='Helvetica'),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = NA, colour = NA), 
    panel.border = element_rect(fill = NA, colour = NA),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )

cowplot::plot_grid(
  cowplot::plot_grid(plot.geneCat, plot.geneLocs, plot.repeats, ncol = 3, align = 'hv', axis = 'tblr', rel_widths = c(.3,.3,.2)),
  plot.repeatsOverview, nrow = 2, align = 'h', axis = 'tb', rel_heights = c(.4,.3)
)


# Figure 3 - Bodymap Expression -------------------------------------------

# Bodymap 2.0 files.
files.bam <- list.files('/mnt/data/ccbc_environment/project/hCRISPR/externalData/bodymap/', full.names = T, pattern = 'chr.*bam$')

# Count overlapping uniquely-mapped reads.
files.bam <- Rsamtools::BamFileList(files.bam)
flag.bam <- Rsamtools::scanBamFlag(isSecondaryAlignment = F, isSupplementaryAlignment = F)
params <- Rsamtools::ScanBamParam(which = data.SRSR.GRanges.hg19, flag = flag.bam)
bodymap.overlap <- tibble::as_tibble(Rsamtools::countBam(files.bam, param = params))

# Convert to GRanges.
bodymap.overlap.GRanges <- makeGRangesFromDataFrame(bodymap.overlap, seqnames.field = 'space', keep.extra.columns = T)

# Add original query information.
x <- IRanges::findOverlaps(bodymap.overlap.GRanges, data.SRSR.GRanges.hg19, minoverlap = 4, select = 'first')

bodymap.overlap.GRanges$chr <- as.character(seqnames(bodymap.overlap.GRanges))
bodymap.overlap.GRanges$CRISPR_Id <- data.SRSR.GRanges.hg19[x,]$CRISPR_Id
bodymap.overlap.GRanges$widthSRSR <- width(bodymap.overlap.GRanges)

# Calculate TPM
plotData.bodymap <- tibble::as_tibble(mcols(bodymap.overlap.GRanges)) %>% 
  dplyr::mutate(tissue = factor(gsub('_', ' ', Hmisc::capitalize(trimws(gsub('\\..*', '', gsub('.*Map\\.', '', file))))))) %>% 
  dplyr::group_by(tissue) %>% 
  dplyr::mutate(
    RPK = records / (widthSRSR / 1e3),
    perMillionScale = (1 / sum(RPK)),
    TPM = (RPK * perMillionScale) * 1e6
  ) %>% 
  dplyr::ungroup()

# All samples should add up to 1 million.
plotData.bodymap %>% dplyr::group_by(tissue) %>% dplyr::summarise(sum(TPM))

options(scipen = 10)
plot1 <- ggplot(plotData.bodymap %>% dplyr::filter(TPM >= 5) %>% dplyr::group_by(tissue) %>% dplyr::summarise(Freq = length(unique(CRISPR_Id))), aes(x = tissue, y = Freq, label = Freq, fill = tissue)) + 
  geom_bar(stat = 'identity', color = 'black', width = .8) +
  scale_y_continuous(limits = c(0, 4000)) + 
  scale_fill_manual(values = unname(hues::iwanthue(16)), guide = F) +
  labs(x = NULL, y = '# SRSR elements\nwith TPM > 5') +
  geom_text(color='black', size = 2.5, vjust = -1) +
  theme(
    legend.position = 'bottom', 
    axis.title.y = element_text(size = 8), 
    text = element_text(size=10, family='Helvetica'),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 1),
    panel.background = element_rect(fill = NA, colour = NA), 
    panel.border = element_rect(fill = NA, colour = NA),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )

# Plot expression.
plot2 <- ggplot(plotData.bodymap %>% dplyr::filter(TPM >= 5), aes(x = tissue, y = TPM)) + 
  ggbeeswarm::geom_quasirandom(aes(color = tissue), dodge.width = .8, alpha = .1) +
  stat_boxplot(geom ='errorbar', position = position_dodge(.8), width = .25, alpha = .33) + 
  geom_boxplot(width = .5, position = position_dodge(.8), outlier.shape = NA) +
  scale_color_manual(values = unname(hues::iwanthue(16)), guide = F) +
  scale_y_continuous(trans = scales::pseudo_log_trans(), limits = c(4, 1000000), breaks = c(10, 25, 50, 100, 250, 500, 1000, 10000, 100000, 1000000)) +
  labs(x = 'BodyMap 2.0 tissues', y = 'Transcripts Per Million (TPM)\n(based on uniquely-mapped reads; log10)') + 
  theme(
    legend.position = 'bottom', 
    text = element_text(size = 8, family='Helvetica', colour = 'black'),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_blank(),
    axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
    panel.background = element_rect(fill = NA, colour = NA)
  )

cowplot::plot_grid(plot1, plot2, align = 'hv', axis = 'tblr', nrow = 2, rel_heights = c(.15, .3))


# Heatmap BodyMap ---------------------------------------------------------

z <- (plotData.bodymap %>% dplyr::filter(TPM >= 100) %>% distinct(CRISPR_Id))$CRISPR_Id
x <- reshape2::dcast(plotData.bodymap %>% dplyr::filter(CRISPR_Id %in% z), CRISPR_Id~tissue, fun.aggregate = sum, value.var = 'TPM')
rownames(x) <- x$CRISPR_Id; x$CRISPR_Id <- NULL

pheatmap::pheatmap(
  x, 
  show_rownames = F,
  color=colorRampPalette(c('green', 'white', 'red'))(50), 
  legend = T, angle_col = 45, cellwidth = 12,
  fontsize = 8, treeheight_col = 25, treeheight_row = 25,
  scale = 'row', clustering_distance_cols = 'euclidean', clustering_distance_rows = 'euclidean',
  clustering_method = 'ward.D2'
)


# Example of Bodymap Expression -------------------------------------------

# Additional info.
data.ENCODE.DNAse <- makeGRangesFromDataFrame(readr::read_delim('/mnt/data/ccbc_environment/project/hCRISPR/externalData/ENCODE/wgEncodeRegDnaseClusteredV3_hg19.bed', delim = '\t', trim_ws = T, col_names = c('chr', 'start', 'end', 'name', 'score', 'expCount',	'expNums',	'expScores')))
data.ENCODE.CpG <- makeGRangesFromDataFrame(readr::read_delim('/mnt/data/ccbc_environment/project/hCRISPR/externalData/ENCODE//cpgIslandExt.hg19.bed', delim = '\t', trim_ws = T, col_names = c('chr', 'start', 'end', 'name', 'score', 'expCount',	'expNums',	'expScores')))
options(ucscChromosomeNames=FALSE)

# Load genome tracks.
track.ideo <- Gviz::IdeogramTrack(genome = 'hg19', chromosome = 'chr9')
track.ideoAxis <- Gviz::GenomeAxisTrack()

pdf(file = 'SRSR.pdf')
lapply((plotData.bodymap %>% dplyr::filter(TPM >= 1000) %>% distinct(CRISPR_Id))$CRISPR_Id[1:2], function(x){
  SRSR <- data.SRSR.GRanges.hg19[data.SRSR.GRanges.hg19$CRISPR_Id == x]
  
  chrom <- as.character(seqnames(SRSR))
  startPos <- start(SRSR)
  endPos <- end(SRSR)
  id <- SRSR$CRISPR_Id
  
  # Get SRSR and expression tracks.
  track.SRSR <- Gviz::SequenceTrack(genome = 'hg19', chromosome = chrom)
  track.SRSR <- Gviz::AnnotationTrack(start=c(startPos), end = c(endPos), chromosome = chrom, id = id, name='SRSR', fill = '#C90035', border = '#C90035', showID = T, fontcolor = 'black', just.group = 'below', cex = .6)
  
  bodymap.files <- list.files('/mnt/data/ccbc_environment/project/hCRISPR/externalData/bodymap/', full.names = T, pattern = 'bam$')
  track.bodymap <- lapply(bodymap.files, function(f){
    tissue <- Hmisc::capitalize(gsub('\\..*', '', gsub('chr_GRCh37.HumanBodyMap.', '', basename(f))))
    track <- Gviz::AlignmentsTrack(f, isPaired = T, chromosome = chrom, name = tissue, ylim = c(0,100), fill = '#FFAAA5')
    return(track)
  })
  
  # track.adrenal <- Gviz::AlignmentsTrack('/mnt/data/ccbc_environment/project/hCRISPR/externalData/bodymap/chr_GRCh37.HumanBodyMap.adrenal.1.bam', isPaired = T, chromosome = chrom, name = 'Adrenal', ylim = c(0,100), fill = '#A8E6CF')
  # track.brain <- Gviz::AlignmentsTrack('/mnt/data/ccbc_environment/project/hCRISPR/externalData/bodymap/chr_GRCh37.HumanBodyMap.brain.1.bam', isPaired = T, chromosome = chrom, name = 'Brain', ylim = c(0,100), fill = '#FFAAA5')
  # track.breast <- Gviz::AlignmentsTrack('/mnt/data/ccbc_environment/project/hCRISPR/externalData/bodymap/chr_GRCh37.HumanBodyMap.breast.1.bam', isPaired = T, chromosome = chrom, name = 'Breast', ylim = c(0,100), fill = '#95C7D9')
  # track.prostate <- Gviz::AlignmentsTrack('/mnt/data/ccbc_environment/project/hCRISPR/externalData/bodymap/chr_GRCh37.HumanBodyMap.prostate.1.bam', isPaired = T, chromosome = chrom, name = 'Prostate', ylim = c(0,100), fill = '#95C7D9')
  # track.lymph <- Gviz::AlignmentsTrack('/mnt/data/ccbc_environment/project/hCRISPR/externalData/bodymap/chr_GRCh37.HumanBodyMap.lymph.1.bam', isPaired = T, chromosome = chrom, name = 'Lymph', ylim = c(0,100), fill = '#DCEDC1')
  
  data.ENCODE.CpG.overlap <- subsetByOverlaps(data.ENCODE.CpG, GRanges(chrom, ranges = IRanges(startPos-5000, endPos+5000)))
  track.CpG <- Gviz::AnnotationTrack(data.ENCODE.CpG.overlap, name = 'CpG', fill = '#00ff80', border = '#00ff80', showID = F)
  
  data.ENCODE.DNAse.overlap <- subsetByOverlaps(data.ENCODE.DNAse, GRanges(chrom, ranges = IRanges(startPos-5000, endPos+5000)))
  track.DNAse <- Gviz::AnnotationTrack(data.ENCODE.DNAse.overlap, name = 'DNAse', fill = '#0080ff', border = '#0080ff',showID = F)
  
  Gviz::plotTracks(c(list(track.ideo, track.ideoAxis, track.SRSR, track.DNAse, track.CpG), track.adrenal), from = startPos-1000, to= endPos+1000, chromosome = chrom, type = 'coverage', sizes =c(.3,.4,.2,.2,.2,rep(.25, 1)), legend = TRUE, showFeatureId = TRUE, just.group = 'above')
  
})

dev.off()


# Figure 4 - Cell-line expression -----------------------------------------

# Cell-line files.
files.bam <- list.files('/mnt/data/ccbc_environment/project/hCRISPR/externalData/Celllines/untrimmedBAM/', full.names = T, pattern = '.*bam$')

# Count overlapping uniquely-mapped reads.
files.bam <- Rsamtools::BamFileList(files.bam)
flag.bam <- Rsamtools::scanBamFlag(isSecondaryAlignment = F, isSupplementaryAlignment = F)
params <- Rsamtools::ScanBamParam(which = dropSeqlevels(data.SRSR.GRanges.hg19, 'chrY', 'coarse'), flag = flag.bam)
celllines.overlap <- tibble::as_tibble(Rsamtools::countBam(files.bam, param = params))

# Convert to GRanges.
celllines.overlap.GRanges <- makeGRangesFromDataFrame(celllines.overlap, seqnames.field = 'space', keep.extra.columns = T)

# Add original query information.
x <- IRanges::findOverlaps(celllines.overlap.GRanges, data.SRSR.GRanges.hg19, minoverlap = 4, select = 'first')

celllines.overlap.GRanges$chr <- as.character(seqnames(celllines.overlap.GRanges))
celllines.overlap.GRanges$CRISPR_Id <- data.SRSR.GRanges.hg19[x,]$CRISPR_Id
celllines.overlap.GRanges$widthSRSR <- width(celllines.overlap.GRanges)

# Calculate TPM
plotData.celllines <- tibble::as_tibble(mcols(celllines.overlap.GRanges)) %>% 
  dplyr::mutate(tissue = factor(gsub('.bam', '', gsub('Nucleus.*Aln', '', gsub('Cell.*Aln', '', gsub('.*RnaSeq', '', celllines.overlap.GRanges$file)))))) %>% 
  dplyr::group_by(tissue) %>% 
  dplyr::mutate(
    RPK = records / (widthSRSR / 1e3),
    perMillionScale = (1 / sum(RPK)),
    TPM = (RPK * perMillionScale) * 1e6
  ) %>% 
  dplyr::ungroup()

levels(plotData.celllines$tissue) <- c('HeLa S3 - Rep #1', 'HeLa S3 - Rep #2', 'Hep G2', 'HUV-EC-C - Rep #3', 'HUV-EC-C - Rep #4', 'K562')

# All samples should add up to 1 million.
plotData.celllines %>% dplyr::group_by(tissue) %>% dplyr::summarise(sum(TPM))

# Total count.
options(scipen = 10)
plot1 <- ggplot(plotData.celllines %>% dplyr::filter(TPM >= 5) %>% dplyr::group_by(tissue) %>% dplyr::summarise(Freq = length(unique(CRISPR_Id))), aes(x = tissue, y = Freq, label = Freq, fill = tissue)) + 
  geom_bar(stat = 'identity', color = 'black', width = .8) +
  scale_y_continuous(limits = c(0, 800)) + 
  scale_fill_manual(values = c('#ffb6b9', '#bbded6', '#fae3d9', '#8ac6d1', '#cd6684', '#6f5a7e'), guide = F) +
  labs(x = 'Cell-lines', y = '# SRSR elements\nwith TPM > 5') +
  geom_text(color='black', size = 2.5, vjust = -1) +
  theme(
    legend.position = 'bottom', 
    axis.title.y = element_text(size = 8), 
    text = element_text(size=10, family='Helvetica'),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 1),
    panel.background = element_rect(fill = NA, colour = NA), 
    panel.border = element_rect(fill = NA, colour = NA),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )

# Plot expression.
plot2 <- ggplot(plotData.celllines %>% dplyr::filter(TPM >= 5), aes(x = tissue, y = TPM)) + 
  ggbeeswarm::geom_quasirandom(aes(color = tissue), dodge.width = .8, alpha = .33) +
  stat_boxplot(geom ='errorbar', position = position_dodge(.8), width = .25, alpha = .33) + 
  geom_boxplot(width = .5, position = position_dodge(.8), outlier.shape = NA) +
  scale_color_manual(values = c('#ffb6b9', '#bbded6', '#fae3d9', '#8ac6d1', '#cd6684', '#6f5a7e'), guide = F) +
  scale_y_continuous(trans = scales::pseudo_log_trans(), limits = c(9, 1000000), breaks = c(10, 25, 50, 100, 250, 500, 1000, 10000, 100000, 1000000)) +
  labs(x = 'Cell-lines', y = 'Transcripts Per Million (TPM)\n(based on uniquely-mapped reads; log10)') + 
  theme(
    legend.position = 'bottom', 
    text = element_text(size = 8, family='Helvetica', colour = 'black'),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_blank(),
    axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
    panel.background = element_rect(fill = NA, colour = NA)
  )

cowplot::plot_grid(plot1, plot2, align = 'hv', axis = 'tblr', nrow = 2, rel_heights = c(.15, .3))

# Overlap with ENCODE TF, Cpg, DNASeI -------------------------------------

#######################
# Downloaded from UCSC Browser (hg38)
#######################

# Transcription factors.
data.ENCODE.TF <- readr::read_delim('/mnt/data/ccbc_environment/project/hCRISPR/externalData/ENCODE/wgEncodeRegTfbsClustered_hg38.bed', delim = '\t', trim_ws = T)
data.ENCODE.TF <- GenomicRanges::makeGRangesFromDataFrame(data.ENCODE.TF, keep.extra.columns = T)
data.ENCODE.TF <- R2CCBC::cleanSeqlevels(data.ENCODE.TF, genome = 'hg38')
data.ENCODE.TF.PerFactor <- GRangesList(split(data.ENCODE.TF, data.ENCODE.TF$name))

# DNAse clusters
data.ENCODE.DNAse <- readr::read_delim('/mnt/data/ccbc_environment/project/hCRISPR/externalData/ENCODE/wgEncodeRegDnaseClustered_hg38.bed', delim = '\t', trim_ws = T)
data.ENCODE.DNAse <- GenomicRanges::makeGRangesFromDataFrame(data.ENCODE.DNAse, keep.extra.columns = T)
data.ENCODE.DNAse <- R2CCBC::cleanSeqlevels(data.ENCODE.DNAse, genome = 'hg38')
data.ENCODE.DNAse$name <- 'DNAse clusters'

# CpG clusters
data.ENCODE.CpG <- readr::read_delim('/mnt/data/ccbc_environment/project/hCRISPR/externalData/ENCODE/cpgIsland_hg38.bed', delim = '\t', trim_ws = T)
data.ENCODE.CpG <- GenomicRanges::makeGRangesFromDataFrame(data.ENCODE.CpG, keep.extra.columns = T)
data.ENCODE.CpG <- R2CCBC::cleanSeqlevels(data.ENCODE.CpG, genome = 'hg38')
data.ENCODE.CpG$name <- 'CpG islands'

# Function to find overlap.
findOverlapData <- function(SRSR, x){
  tibble::as_tibble(data.frame(
    overlapSRSR = length(unique(queryHits(findOverlaps(SRSR, unique(x), minoverlap = 5)))),
    factor = unique(x$name),
    totalPeaksInFactor = length(unique(x))
  )) %>% dplyr::mutate(
    overlapSRSR.relative = overlapSRSR / length(SRSR)
  )
}

# Overlap TF.
data.ENCODE.TF.PerFactor.overlap <- do.call(rbind, lapply(data.ENCODE.TF.PerFactor, function(x){
  findOverlapData(data.SRSR.GRanges, x)
}))

# Overlap DNAse
data.ENCODE.DNAse.overlap <- findOverlapData(data.SRSR.GRanges, data.ENCODE.DNAse)

# Overlap CpG
data.ENCODE.CpG.overlap <- findOverlapData(data.SRSR.GRanges, data.ENCODE.CpG)

# Plot TF.
plotData <- data.ENCODE.TF.PerFactor.overlap %>% dplyr::arrange(-overlapSRSR.relative) %>% 
  head(15) %>%
  reshape2::melt(., id.vars = c('factor', 'totalPeaksInFactor')) %>% 
  dplyr::filter(grepl('relative', variable)) %>% 
  dplyr::mutate(identifier = sprintf('%s\n(%s peaks)', factor, totalPeaksInFactor))

plot1 <- ggplot(plotData, aes(x = reorder(identifier, value), y = value, fill = variable, label = paste0(round(value, 3) * 100, '%'))) + 
  geom_bar(stat = 'identity', position = 'dodge', fill = '#005288', color = 'black') +
  scale_y_continuous(limits = c(0, .06), labels = scales::percent) +
  geom_text(color='black', size = 2.5, hjust = -.15) +
  labs(x = 'Transcription factors\n(Top 15)', y = 'Rel. overlap of SRSR elements\nwith genomic binding sites') + 
  theme(
    legend.position = 'bottom', 
    axis.title.y = element_text(size = 8), 
    text = element_text(size=10, family='Helvetica'),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.y = element_line(colour = 'grey80', linetype = 'dotted'),
    panel.background = element_rect(fill = NA, colour = NA), 
    panel.border = element_rect(fill = NA, colour = NA),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  ) + coord_flip()

# Plot CpG and DNAse.
plotData <- rbind(data.ENCODE.DNAse.overlap, data.ENCODE.CpG.overlap) %>% 
  reshape2::melt(., id.vars = c('factor', 'totalPeaksInFactor')) %>% 
  dplyr::filter(grepl('relative', variable)) %>% 
  dplyr::mutate(identifier = sprintf('%s\n(%s peaks)', factor, totalPeaksInFactor))

plot2 <- ggplot(plotData, aes(x = reorder(identifier, -value), y = value, fill = variable, label = paste0(round(value, 3) * 100, '%'))) + 
  geom_bar(stat = 'identity', position = 'dodge', fill = '#ED7A8B', color = 'black') +
  scale_y_continuous(limits = c(0, .25), labels = scales::percent) +
  geom_text(color='black', size = 2.5, vjust = -1) +
  labs(x = NULL, y = 'Rel. overlap of SRSR elements\nwith genomic locations') + 
  theme(
    legend.position = 'bottom', 
    axis.title.y = element_text(size = 8), 
    text = element_text(size=10, family='Helvetica'),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = NA, colour = NA), 
    panel.border = element_rect(fill = NA, colour = NA),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )

# Combine plots.
cowplot::plot_grid(plot1, plot2, align = 'hv', axis = 'tblr', ncol = 2, rel_widths = c(.2,.15))



# Overlap with DASHR v2.0, piRBase & RNACentral  --------------------------

# Import DASHR v2.0 data (hg38).
dashr2 <- readr::read_delim('/mnt/data/ccbc_environment/project/hCRISPR/externalData/Dashr2/dashr.v2.annotation.hg38.gff', delim = '\t', col_names = F)
dashr2 <- dashr2 %>% dplyr::select(X1, X2, X3, X4, X5, X7, X9)
dashr2 <- R2CCBC::cleanSeqlevels(makeGRangesFromDataFrame(dashr2, seqnames.field = 'X1', start.field = 'X4', end.field = 'X5', strand.field = 'X7', keep.extra.columns = T), genome = 'hg38')
colnames(mcols(dashr2)) <- c('Source', 'Type', 'ID')
dashr2$id <- sprintf('%s (%s)', gsub('ID\\=', '', gsub(';.*', '', dashr2$ID)), dashr2$Type)

# Import piRBase data (hg38).
piRBase <- rtracklayer::import.bed('/mnt/data/ccbc_environment/project/hCRISPR/externalData/piRBase/hsa.bed.gz')
piRBase <- R2CCBC::cleanSeqlevels(piRBase, genome = 'hg38')

# Import RNAcentral (hg38; v15)
rnacentral <- rtracklayer::import.bed('/mnt/data/ccbc_environment/project/hCRISPR/externalData/RNACentral/homo_sapiens.GRCh38.bed')
rnacentral <- R2CCBC::cleanSeqlevels(rnacentral, genome = 'hg38')
rnacentral$name <- sprintf('%s (%s)', rnacentral$name, rnacentral$NA.1)

# Find overlap.
overlap.dashr2 <- data.frame(mcols(findOverlapSRSR(data.SRSR.GRanges, dashr2))) %>% dplyr::group_by(CRISPR_Id) %>% dplyr::summarise(dashr2 = paste0(unique(id), collapse = ', '), type = paste0(unique(Type), collapse = ', ')) %>% tidyr::separate_rows(type, sep = ", ", convert = FALSE)
overlap.piRBase <- data.frame(mcols(findOverlapSRSR(data.SRSR.GRanges, piRBase))) %>% dplyr::group_by(CRISPR_Id) %>% dplyr::summarise(pirBase = paste0(unique(name), collapse = ', '), type = 'piRNA')
overlap.rnacentral <- data.frame(mcols(findOverlapSRSR(data.SRSR.GRanges, rnacentral))) %>% dplyr::group_by(CRISPR_Id) %>% dplyr::summarise(rnacentral = paste0(unique(name), collapse = ', '), type = paste0(unique(NA.1), collapse = ', ')) %>% tidyr::separate_rows(type, sep = ", ", convert = FALSE)

plotData <- rbind(
  overlap.dashr2 %>% dplyr::select(CRISPR_Id, type),
  overlap.piRBase %>% dplyr::select(CRISPR_Id, type),
  overlap.rnacentral%>% dplyr::select(CRISPR_Id, type)
) %>% 
  dplyr::mutate(type = trimws(gsub('_', '', type))) %>% 
  dplyr::mutate(type = dplyr::if_else(grepl('misc', type), 'Misc. RNA', type)) %>% 
  dplyr::mutate(type = dplyr::if_else(grepl('SRPRNA', type), 'SRP-RNA', type)) %>% 
  dplyr::distinct(CRISPR_Id, type) %>% 
  dplyr::group_by(type) %>% 
  dplyr::mutate(nType = dplyr::n_distinct(CRISPR_Id)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(type = dplyr::if_else(nType <= 25, 'Misc. RNA', type)) %>% 
  dplyr::group_by(CRISPR_Id) %>% 
  dplyr::summarise(type = type %>% unique %>% sort %>% paste(collapse = " & ")) %>% 
  dplyr::group_by(type) %>% dplyr::summarise(count = dplyr::n_distinct(CRISPR_Id)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(type = dplyr::if_else(count < 100, 'Multiple RNA classes (≥2)\nwith <100 overlapping\nSRSR per combination', type)) %>% 
  dplyr::group_by(type) %>% dplyr::summarise(count = sum(count))

plot.smallRNA <-
  ggplot(plotData, aes(x = reorder(type, -count), y = count, fill = type, label = count)) + 
  geom_bar(stat = 'identity', position = 'dodge', fill = '#ED7A8B', color = 'black') +
  scale_y_continuous(limits = c(0, 7500)) +
  geom_text(color='black', size = 2.5, vjust = -1) +
  labs(x = 'RNA Classes', y = '# overlapping SRSR elements') + 
  theme(
    legend.position = 'bottom', 
    axis.title.y = element_text(size = 8), 
    text = element_text(size=10, family='Helvetica'),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = NA, colour = NA), 
    panel.border = element_rect(fill = NA, colour = NA),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )

# Determine width of overlap.
x <- rnacentral[grepl('piRNA', rnacentral$NA.1),]
mcols(x) <- NULL

x <- unlist(GRangesList(
  x,
  IRanges::subsetByOverlaps(piRBase, data.SRSR.GRanges, minoverlap = 5),
  dashr2[grepl('piRNA', dashr2$Type),]
))

overlapping.index = findOverlaps(query = data.SRSR.GRanges, subject = x, minoverlap = 5)
summary(width(data.SRSR.GRanges[queryHits(overlapping.index)]) / width(x[subjectHits(overlapping.index)]))


# Import EPCATS -----------------------------------------------------------

epcatSites <- rtracklayer::import.bed('/mnt/data/ccbc_environment/project/hCRISPR/externalData/150323_EPCAT_def_with_validated_hg19.bed')
epcatSites <- R2CCBC::cleanSeqlevels(epcatSites)


# Import U133 -------------------------------------------------------------

U133.data <- rtracklayer::import.bed('/mnt/data/ccbc_environment/project/hCRISPR/externalData/U133/HG-U133_Plus_2.hg19.bed')
u133.probes <- R2CCBC::cleanSeqlevels(U133.data$`Affymetrix HG-U133_Plus_2 Probe`)
u133.regions <- R2CCBC::cleanSeqlevels(U133.data$`Affymetrix HG-U133_Plus_2 Consensus/Exemplar`)


# Perform overlap ---------------------------------------------------------

# Overlap per SRSR.
overlap.humanRepeats <- data.frame(mcols(findOverlapSRSR(data.SRSR.GRanges, humanRepeats))) %>% dplyr::group_by(CRISPR_Id) %>% dplyr::summarise(repeatFamily = paste0(unique(repeatFamily), collapse = ', '))
overlap.ENCODE.CpG <- data.frame(mcols(findOverlapSRSR(data.SRSR.GRanges, data.ENCODE.CpG))) %>% dplyr::group_by(CRISPR_Id) %>% dplyr::summarise(CpG = 'Yes')
overlap.ENCODE.DNAse <- data.frame(mcols(findOverlapSRSR(data.SRSR.GRanges, data.ENCODE.DNAse))) %>% dplyr::group_by(CRISPR_Id) %>% dplyr::summarise(DNAse = 'Yes')
overlap.ENCODE.TF <- data.frame(mcols(findOverlapSRSR(data.SRSR.GRanges, data.ENCODE.TF))) %>% dplyr::group_by(CRISPR_Id) %>% dplyr::summarise(TF = paste0(unique(name), collapse = ', '))
overlap.HumanGenes <- data.frame(mcols(findOverlapSRSR(data.SRSR.GRanges, geneAnnotation))) %>% dplyr::group_by(CRISPR_Id) %>% dplyr::summarise(Gene = paste0(unique(gene_name), collapse = ', '), GeneType = paste0(unique(gene_type), collapse = ', '))
overlap.EPCATS <- data.frame(mcols(findOverlapSRSR(data.SRSR.GRanges.hg19, epcatSites))) %>% dplyr::group_by(CRISPR_Id) %>% dplyr::summarise(epcats = paste0(unique(name), collapse = ', '))
overlap.U133.probes <- data.frame(mcols(findOverlapSRSR(data.SRSR.GRanges.hg19, u133.probes))) %>% dplyr::group_by(CRISPR_Id) %>% dplyr::summarise(u133.probes = paste0(unique(name), collapse = ', '))
overlap.U133.regions <- data.frame(mcols(findOverlapSRSR(data.SRSR.GRanges.hg19, u133.regions))) %>% dplyr::group_by(CRISPR_Id) %>% dplyr::summarise(u133.regions = paste0(unique(name), collapse = ', '))
overlap.bodymap <- plotData.bodymap %>% dplyr::group_by(CRISPR_Id, tissue) %>%  dplyr::summarise(TPM = sum(TPM)) %>% reshape2::dcast(CRISPR_Id ~ tissue)
overlap.celllines <- plotData.celllines %>% dplyr::group_by(CRISPR_Id, tissue) %>%  dplyr::summarise(TPM = sum(TPM)) %>% reshape2::dcast(CRISPR_Id ~ tissue)
overlap.dashr2
overlap.piRBase


# Output all data ---------------------------------------------------------

data.SRSR.tibble <- tibble::as_tibble(data.frame(data.SRSR.GRanges))
data.SRSR.tibble <- data.SRSR.tibble %>% 
  dplyr::left_join(overlap.humanRepeats) %>% 
  dplyr::left_join(overlap.ENCODE.CpG) %>% 
  dplyr::left_join(overlap.ENCODE.DNAse) %>% 
  dplyr::left_join(overlap.ENCODE.TF) %>% 
  dplyr::left_join(overlap.HumanGenes) %>% 
  dplyr::left_join(overlap.EPCATS) %>% 
  dplyr::left_join(overlap.U133.probes) %>% 
  dplyr::left_join(overlap.U133.regions) %>% 
  dplyr::left_join(overlap.bodymap) %>% 
  dplyr::left_join(overlap.celllines) %>% 
  dplyr::left_join(overlap.dashr2) %>% 
  dplyr::left_join(overlap.piRBase) %>% 
  dplyr::left_join(overlap.rnacentral)

write.table(data.SRSR.tibble, file = 'SupplTable1.txt', append = F, quote = F, row.names = F, sep = '\t')
