# Author:                      Job van Riet
# Date of  creation:           05-02-2020
# Function:                    Generating suppl. figures of SRSR manuscript.


# Packages and general functions ------------------------------------------

library(plyr)
library(dplyr)
library(ggplot2)


# Import data -------------------------------------------------------------

# Import hg38 CRISPRCASFinder results.
data.SRSR <- readr::read_delim('/mnt/data/ccbc_environment/project/hCRISPR/CRISPRCasFinder/hg38_2019_11_29/TSV/Crisprs_REPORT.tsv', delim = '\t', trim_ws = T)

# Convert to GRanges.
data.SRSR.GRanges <- GenomicRanges::GRanges(seqnames = data.SRSR$Sequence, ranges = IRanges(data.SRSR$CRISPR_Start, data.SRSR$CRISPR_End))
GenomicRanges::mcols(data.SRSR.GRanges) <- data.SRSR
seqlevelsStyle(data.SRSR.GRanges) = "UCSC"

# Remove haplo-chromosomes, if any.
data.SRSR.GRanges <- R2CCBC::cleanSeqlevels(GenomicRanges::sort(data.SRSR.GRanges), genome ='hg38')

# Make Tibble for plotting.
data.SRSR.tibble <- tibble::as_tibble(data.SRSR.GRanges)

# Plot genomic characterics -----------------------------------------------

# Number of SRSR per chromosome.
plotA <- ggplot(data.SRSR.tibble, aes(x = seqnames, fill = seqnames)) + 
  geom_bar(position = 'dodge', width = .75, color = 'black') +
  geom_text(stat= 'count', aes(label = ..count..), position = position_dodge(width = .75), vjust=-.75, size = 2.5) +
  labs(y = '# SRSR') + 
  scale_fill_manual(values = unname(hues::iwanthue(24)), guide = F) +
  theme(
    legend.position = 'right',
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 8),
    text=element_text(size=10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'grey20', linetype = 'dotted'),
    panel.grid.minor.y = element_line(colour = 'grey50', linetype = 'dotted'),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )

# Length of SRSR
plotB <- ggplot(data.SRSR.tibble, aes(x = seqnames, y = CRISPR_Length)) + 
  ggbeeswarm::geom_quasirandom(aes(color = seqnames), dodge.width = .8, alpha = .1) +
  stat_boxplot(geom ='errorbar', position = position_dodge(.8), width = .25, alpha = .33) + 
  geom_boxplot(width = .5, position = position_dodge(.8), outlier.shape = NA) +
  scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 10, 50, 100, 250, 500, 1000, 2000), limits = c(40, 2000)) +
  labs(x = 'Chromosomes', y = 'SRSR element length (bp; log10)') + 
  scale_color_manual(values = unname(hues::iwanthue(24)), guide = F) +
  theme(
    legend.position = 'right',
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 8),
    text=element_text(size=10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'grey20', linetype = 'dotted'),
    panel.grid.minor.y = element_line(colour = 'grey50', linetype = 'dotted'),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )

# Length of Direct Repeat Sequence.
plotC <- ggplot(data.SRSR.tibble, aes(x = seqnames, y = Repeat_Length)) + 
  ggbeeswarm::geom_quasirandom(aes(color = seqnames), dodge.width = .8, alpha = .1) +
  stat_boxplot(geom ='errorbar', position = position_dodge(.8), width = .25, alpha = .33) + 
  geom_boxplot(width = .5, position = position_dodge(.8), outlier.shape = NA) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60), limits = c(18, 60)) +
  labs(x = 'Chromosomes', y = 'DR length of SRSR (bp)') + 
  scale_color_manual(values = unname(hues::iwanthue(24)), guide = F) +
  theme(
    legend.position = 'right',
    axis.text.x = element_text(angle =45, hjust = 1),
    axis.title.y = element_text(size = 8),
    text=element_text(size=10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'grey20', linetype = 'dotted'),
    panel.grid.minor.y = element_line(colour = 'grey50', linetype = 'dotted'),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )

cowplot::plot_grid(plotA, plotB, plotC, nrow = 3, align = 'hv', axis = 'tblr', rel_heights = c(.33, .5, .5), labels = letters[1:3])


# Trend to chromosome length ----------------------------------------------

x <- data.frame(chrLength = GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
x$seqnames <- rownames(x)

corData <- data.SRSR.tibble %>% dplyr::group_by(seqnames) %>% dplyr::tally()
corData <- corData %>% dplyr::inner_join(x)

plotD <- ggplot(corData, aes(x = n, y = chrLength, label = seqnames)) + 
  scale_y_continuous(limits = c(0, 300000000)) +
  geom_point() + 
  ggpmisc::stat_poly_eq(formula = x~y,  aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  parse = TRUE) + 
  geom_smooth(method="lm", se = T, alpha = .1) +
  ggrepel::geom_text_repel(aes(color = seqnames), size = 3) +
  labs(x = '# SRSR elements', y = 'Chromosome length (bp)') +
  scale_color_manual(values = unname(hues::iwanthue(24)), guide = F) +
  theme(
    legend.position = 'right',
    axis.title.y = element_text(size = 8),
    text = element_text(size=10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'grey20', linetype = 'dotted'),
    panel.grid.minor.y = element_line(colour = 'grey50', linetype = 'dotted'),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = 'bold'),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'black', fill = c('darkgreen', 'tomato'))
  )


# Re-visualize U133 markers -----------------------------------------------

data.TB <- readxl::read_excel('/mnt/data/ccbc_environment/project/hCRISPR/externalData/Suppl Table 8 Disease Related Probes U133.xlsx', sheet = 1, trim_ws = T, col_names = T) %>% reshape2::melt() %>% dplyr::mutate(variable = gsub('\\..*', '', variable), Disease = 'Tuberculosis')
data.Asthma <- readxl::read_excel('/mnt/data/ccbc_environment/project/hCRISPR/externalData//Suppl Table 8 Disease Related Probes U133.xlsx', sheet = 2, trim_ws = T, col_names = T) %>% reshape2::melt() %>% dplyr::mutate(variable = gsub('\\..*', '', variable), Disease = 'Asthma')
data.CBLL <- readxl::read_excel('/mnt/data/ccbc_environment/project/hCRISPR/externalData//Suppl Table 8 Disease Related Probes U133.xlsx', sheet = 3, trim_ws = T, col_names = T) %>% reshape2::melt() %>% dplyr::mutate(variable = gsub('\\..*', '', variable), Disease = 'CBLL')
data.Franciscella <- readxl::read_excel('/mnt/data/ccbc_environment/project/hCRISPR/externalData//Suppl Table 8 Disease Related Probes U133.xlsx', sheet = 4, trim_ws = T, col_names = T) %>% reshape2::melt() %>% dplyr::mutate(variable = gsub('\\.\\..*', '', variable), Disease = 'Franciscella Infection')
data.HCC <- readxl::read_excel('/mnt/data/ccbc_environment/project/hCRISPR/externalData//Suppl Table 8 Disease Related Probes U133.xlsx', sheet = 5, trim_ws = T, col_names = T) %>% reshape2::melt() %>% dplyr::mutate(variable = gsub('\\.\\..*', '', variable), Disease = 'Hepatocellular Carcinoma')
data.NDPIT <- readxl::read_excel('/mnt/data/ccbc_environment/project/hCRISPR/externalData//Suppl Table 8 Disease Related Probes U133.xlsx', sheet = 6, trim_ws = T, col_names = T) %>% reshape2::melt() %>% dplyr::mutate(variable = gsub('\\.\\..*', '', variable), Disease = 'NDPIT vs. CDPIT')
data.PCa <- readxl::read_excel('/mnt/data/ccbc_environment/project/hCRISPR/externalData//Suppl Table 8 Disease Related Probes U133.xlsx', sheet = 7, trim_ws = T, col_names = T) %>% reshape2::melt() %>% dplyr::mutate(variable = gsub('\\.\\..*', '', variable), Disease = 'Prostate cancer')
data.Teratozoospermia <- readxl::read_excel('/mnt/data/ccbc_environment/project/hCRISPR/externalData//Suppl Table 8 Disease Related Probes U133.xlsx', sheet = 8, trim_ws = T, col_names = T) %>% reshape2::melt() %>% dplyr::mutate(variable = gsub('\\.\\..*', '', variable), Disease = 'Teratozoospermia')

cowplot::plot_grid(
  makeBox(data.TB, c(5, 10)),
  makeBox(data.Asthma, c(20, 37.5)),
  makeBox(data.CBLL, c(-2, 2)),
  makeBox(data.Franciscella, c(0, 8.5)),
  makeBox(data.HCC, c(2.5, 7.5)),
  makeBox(data.NDPIT, c(250, 1250)),
  makeBox(data.PCa, c(0, 7500)),
  makeBox(data.Teratozoospermia, c(0, 150)),
  ncol = 4, align = 'hv', axis = 'tblr'
)

makeBox <- function(data, ylim){
  ggplot(data, aes(variable, value, fill = variable)) + 
    scale_y_continuous(limits = ylim) +
    stat_boxplot(geom ='errorbar', width = .25, alpha = .33) +
    geom_boxplot(width = .5, notch = F, outlier.shape = NA, alpha = .33) +
    ggbeeswarm::geom_quasirandom(dodge.width = .5, cex = 1.5, alpha = 1) +
    scale_fill_manual( values = c('darkblue', '#ff6145', 'orange'), guide = guide_legend(title = 'Disease state', title.position = 'top', title.hjust = 0.5, nrow = 2, keywidth = 1, keyheight = 1)) +
    ggsignif::stat_signif(comparisons=combn(unique(data$variable), m = 2, simplify = F), show.legend = F, test = 'wilcox.test', map_signif_level = T, step_increase = .04, color = 'black', tip_length = .01) +
    # Axis titles.
    labs(x = NULL, y = 'Transcript Expression Intensity', subtitle = sprintf('U133 Probe: %s', unique(data$Probe))) + 
    # Theme
    theme(legend.position = 'bottom', 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          text=element_text(size=10, family='Helvetica'),
          panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
          panel.grid.major.x = element_blank(),
          panel.background = element_rect(fill = 'white', colour = NA),
          panel.border = element_rect(fill = NA, colour = 'grey20')
    )
  
}
