# Author:                      Job van Riet
# Date of  creation:           19-06-2020
# Function:                    Analysis short-read sequencing.

# Packages and general functions ------------------------------------------

library(plyr)
library(dplyr)
library(ggplot2)


# Count the reads per SRSR ------------------------------------------------

# for x in *bam; do echo "samtools view -F 256 ${x} | grep -E 'chr' > ${x}.alignedReads"; done


# Load metadata -----------------------------------------------------------

metadata <- readr::read_delim("/mnt/data/ccbc_environment/project/hCRISPR/externalData/GSE80400/metadata.txt", delim = '\t')

# Overview reads per BAM --------------------------------------------------

data.flagstat <- R2CCBC::readFlagstat.sambamba(list.files('/mnt/data/ccbc_environment/project/hCRISPR/externalData/AlignedShortReads/', pattern = 'flagstat', full.names = T))
data.flagstat <- data.flagstat %>% tibble::rownames_to_column() %>% dplyr::select(Run = rowname, readsPassedQC, readsMapped) %>% dplyr::mutate(Run = gsub('\\.bam.*', '', Run), relReads = (readsMapped / readsPassedQC) * 100)
data.flagstat <- data.flagstat %>% dplyr::inner_join(metadata) %>% dplyr::filter(source_name == 'prostate', storage != 'FFPE')

ggplot(data.flagstat, aes(x = Run, y = readsMapped, label = round(relReads, 2), fill = status)) +
  geom_bar(stat = 'identity', color = 'black') +
  geom_text() +
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


# Import SRSR reads per BAM -----------------------------------------------

importReadsSRSR <- function(x){
  z <- readr::read_delim(x, delim = '\t', col_names = c(1,2,'SRSR', 4:20))
  z <- z %>% dplyr::group_by(SRSR) %>% dplyr::summarise(n = length(SRSR))
  colnames(z) <- c('SRSR', 'nReads')
  z$Run <- unique(gsub('\\.bam.*', '', basename(x)))
  
  return(z)  
}


reads.SRSR <- lapply(list.files('/mnt/data/ccbc_environment/project/hCRISPR/externalData/AlignedShortReads/', pattern = 'SRR.*alignedReads', full.names = T), importReadsSRSR)
reads.SRSR <- do.call(rbind, reads.SRSR)

reads.SRSR <- reads.SRSR %>% dplyr::inner_join(metadata) %>% dplyr::filter(source_name == 'prostate', storage != 'FFPE')
reads.SRSR <- reads.SRSR %>% dplyr::inner_join(data.flagstat %>% dplyr::select(Run, readsMapped))
reads.SRSR <- reads.SRSR %>% dplyr::group_by(SRSR, status) %>% dplyr::summarise(meanReads = sum(nReads), meanTotal = sum(readsMapped))

reads.SRSR <- reads.SRSR %>% dplyr::group_by(SRSR) %>% dplyr::summarise(
  p = chisq.test( matrix(c(meanReads[status == 'normal'], meanReads[status == 'tumor'], meanTotal[status == 'normal'], meanTotal[status == 'tumor']), ncol = 2))$p.value,
  meanReads.Normal = meanReads[status == 'normal'],
  meanReads.Tumor = meanReads[status == 'tumor']) %>% dplyr::ungroup()

reads.SRSR$padj <- p.adjust(reads.SRSR$p)
reads.SRSR <- reads.SRSR %>% dplyr::mutate(padj = ifelse(p == 0, NA, padj))


x <- readr::read_delim('R/srsr.txt', delim = '\t')
x %>% dplyr::left_join(reads.SRSR, by = c('SRSR ID' = 'SRSR')) %>% write.table(., file = 'asd.txt', quote = F, row.names = F, sep = '\t')

