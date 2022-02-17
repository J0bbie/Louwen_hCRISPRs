# Author:                      Job van Riet
# Date of  creation:           15-11-2021
# Function:                    Retrieve the spacers predicted by CRISPRCasFinder.

# Packages and general functions ---

library(plyr)
library(dplyr)


# Import spacers from CRISPRCasFinder ---

readCRISPRCasFinderGFF <- function(file){
  x <- readr::read_tsv(file, skip_empty_rows = T, skip = 1, col_names = c('chr', 'tool', 'type', 'start', 'end', 'x', 'y', 'z', 'info'), trim_ws = T)
  
  x <- x %>% 
    dplyr::filter(type == 'CRISPRspacer') %>% 
    tidyr::separate(.data$info, sep = ';', into = c('sequence', 'name', 'parent', 'ID')) %>% 
    dplyr::mutate(
      sequence = gsub('sequence=', '', sequence),
      parent = gsub('Parent=', '', parent)
    ) %>% 
    tidyr::separate(.data$parent, sep = '_', into = c('chrParent', 'startParent', 'endParent')) %>% 
    dplyr::mutate(
      SRSR = sprintf('%s;%s-%s', chrParent, startParent, endParent)
    ) %>% 
    dplyr::select(SRSR, sequence)
  
  return(x)
}

data.Spacers <- pbapply::pblapply(list.files(path = 'Misc/', pattern = '.gff', full.names = T), readCRISPRCasFinderGFF, cl = 5)
data.Spacers <- dplyr::bind_rows(data.Spacers) %>% 
  dplyr::group_by(SRSR) %>% 
  dplyr::summarise(spacers = paste(sequence, collapse = ',')) %>% 
  dplyr::ungroup()

write.table(data.Spacers, file = 'asd.txt', append = F, quote = F, sep = '\t', row.names = F)