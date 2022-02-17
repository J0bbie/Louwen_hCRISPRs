# Author:                      Job van Riet
# Date of  creation:           15-11-2021
# Function:                    Determine statistical significance of overlap of hCRISPRs with genomic regions.

# Packages and general functions ----

set.seed(123)
library(GenomicRanges)
library(regioneR)

# Import hCRISPRs ----

hCRISPR <- readxl::read_xlsx('/mnt/onco0002/users/jvanriet/hCRISPR/externalData/Supplementary Table 1.xlsx', skip = 3)
hCRISPR <- GRanges(hCRISPR$GRCh38)

# Import genomic regions. ----

# Genomic regions from UCSC Browser (hg38) ----
## Transcription factors.
data.ENCODE.TF <- readr::read_delim('/mnt/onco0002/users/jvanriet/hCRISPR/externalData/ENCODE/wgEncodeRegTfbsClustered_hg38.bed', delim = '\t', trim_ws = T)
data.ENCODE.TF <- GenomicRanges::makeGRangesFromDataFrame(data.ENCODE.TF, keep.extra.columns = T)
data.ENCODE.TF <- keepStandardChromosomes(data.ENCODE.TF, pruning.mode = 'coarse')
data.ENCODE.TF.PerFactor <- GRangesList(split(data.ENCODE.TF, data.ENCODE.TF$name))

## DNAse clusters
data.ENCODE.DNAse <- readr::read_delim('/mnt/onco0002/users/jvanriet/hCRISPR/externalData/ENCODE/wgEncodeRegDnaseClustered_hg38.bed', delim = '\t', trim_ws = T)
data.ENCODE.DNAse <- GenomicRanges::makeGRangesFromDataFrame(data.ENCODE.DNAse, keep.extra.columns = T)
data.ENCODE.DNAse <- keepStandardChromosomes(data.ENCODE.DNAse, pruning.mode = 'coarse')
data.ENCODE.DNAse$name <- 'DNAse clusters'

## CpG clusters
data.ENCODE.CpG <- readr::read_delim('/mnt/onco0002/users/jvanriet/hCRISPR/externalData/ENCODE/cpgIsland_hg38.bed', delim = '\t', trim_ws = T)
data.ENCODE.CpG <- GenomicRanges::makeGRangesFromDataFrame(data.ENCODE.CpG, keep.extra.columns = T)
data.ENCODE.CpG <- keepStandardChromosomes(data.ENCODE.CpG, pruning.mode = 'coarse')
data.ENCODE.CpG$name <- 'CpG islands'


# Determine overlap with hCRISPRs ----

overlapTests <- list()
overlapTests$CpG <- regioneR::overlapPermTest(A = hCRISPR, B = data.ENCODE.CpG, ntimes = 5000, alternative = 'auto', mc.set.seed = FALSE, force.parallel = TRUE)
overlapTests$DNAse <- regioneR::overlapPermTest(A = hCRISPR, B = data.ENCODE.DNAse, ntimes = 5000, alternative = 'auto', mc.set.seed = FALSE, force.parallel = TRUE)
overlapTests$POLR2A <- regioneR::overlapPermTest(A = hCRISPR, B = data.ENCODE.TF.PerFactor$POLR2A, ntimes = 5000, alternative = 'auto', mc.set.seed = FALSE, force.parallel = TRUE)
overlapTests$CTCF <- regioneR::overlapPermTest(A = hCRISPR, B = data.ENCODE.TF.PerFactor$CTCF, ntimes = 5000, alternative = 'auto', mc.set.seed = FALSE, force.parallel = TRUE)
