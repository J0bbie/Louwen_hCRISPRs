# Author:                      Job van Riet
# Date of  creation:           10-02-2022
# Function:                    Visualize the qPCR and assay data.

# Packages and general functions ----

library(dplyr)
library(ggplot2)


# WTS data ----

data.RNASeq <- readxl::read_xlsx('Misc./Suppl. Table X - qPCR.xlsx') %>% 
  dplyr::mutate(Type = factor(Type, levels = c('NAP', 'Tumor')))

data.RNASeq %>% 
  ggplot2::ggplot(., ggplot2::aes(x = Type, y = Counts, fill = SampleID, group = Type)) +
  ggplot2::geom_segment(aes(xend = Type, yend = 0)) +
  ggplot2::geom_point(shape = 21, size = 2) +
  ggplot2::scale_y_continuous(breaks = c(0, 250, 500, 750, 1000), limits = c(0, 1000)) + 
  ggplot2::scale_fill_manual(values = c('#FF1F5B', '#00CD6D', '#009ADE', '#FFC51E'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  ggplot2::labs(x = 'Tissue', y = 'Normalized counts') +
  ggplot2::theme(
    legend.position = 'bottom', 
    text = element_text(size = 8, family='Helvetica', colour = 'black'),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_blank(),
    axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
    panel.background = element_rect(fill = NA, colour = NA)
  ) +
  facet_wrap(SampleID ~ Target, ncol = 2)


# Raw qPCR ----

data.qPCR <- readxl::read_xlsx('Misc./Suppl. Table X - qPCR.xlsx', sheet = 2) %>% 
  dplyr::mutate(Type = factor(Type, levels = c('NAP', 'Tumor'))) %>% 
  dplyr::filter(Target != 'Spike-in',  !is.na(Quantity))

## Copies ----

data.qPCR %>% 
  ggplot2::ggplot(., ggplot2::aes(x = Type, y = Quantity, fill = patientID, group = Type)) +
  ggplot2::geom_segment(aes(xend = Type, yend = 0)) +
  ggplot2::geom_point(shape = 21, size = 2) +
  ggplot2::scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), limits = c(0, 10000), trans = scales::pseudo_log_trans()) +
  ggplot2::scale_fill_manual(values = c('#FF1F5B', '#00CD6D', '#009ADE', '#FFC51E'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  ggplot2::labs(x = 'Tissue', y = 'Quantity') +
  ggplot2::theme(
    legend.position = 'bottom', 
    text = element_text(size = 8, family='Helvetica', colour = 'black'),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_blank(),
    axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
    panel.background = element_rect(fill = NA, colour = NA)
  ) +
  ggplot2::facet_wrap(patientID~Target, ncol = 2)

# Control ----

data.qPCR <- readxl::read_xlsx('Misc./Suppl. Table X - qPCR.xlsx', sheet = 2)

data.qPCR %>% 
  ggplot2::ggplot(., ggplot2::aes(x = patientID, y = Quantity, fill = Type)) +
  ggplot2::geom_bar(stat = 'identity', position = ggplot2::position_dodge2(), col = 'black', lwd = .25) +
  ggplot2::geom_errorbar(aes(ymin = ifelse(Quantity - `Quantity SD` < 0, 0, Quantity - `Quantity SD`), ymax = Quantity + `Quantity SD`), position = 'dodge', size = .5, width = .5) +
  ggplot2::scale_fill_manual( values = c(NAP = '#52639E', Tumor = '#EE6D7A'), guide = guide_legend(title = 'Tissue', title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 1, keyheight = 1)) +
  ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 8)) +
  ggplot2::labs(x = NULL, y = 'No. of copies (Quantity)') +
  ggplot2::facet_wrap(.~Target, nrow = 1, scales = 'free_y') +
  ggplot2::theme(
    legend.position = 'bottom',
    legend.direction = 'horizontal',
    text = ggplot2::element_text(size = 9, family = 'Helvetica', face = 'bold'),
    axis.text.x = ggtext::element_markdown(),
    axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
    axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
    strip.text = ggtext::element_textbox_simple(width = NULL, halign = .5),
    panel.grid.major.x = ggplot2::element_line(colour = 'grey90', linetype = 'dotted'),
    panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
    panel.border = ggplot2::element_rect(fill = NA, colour = NA),
    strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.key = ggplot2::element_blank()
  )

# Assay linearity ----

data.assay <- readxl::read_xlsx('Misc./Suppl. Table X - qPCR.xlsx', sheet = 4) %>% 
  dplyr::filter(!is.na(Quantity_Calculated))

data.assay %>% 
  ggplot2::ggplot(., ggplot2::aes(x = Quantity_Expected, y = Quantity_Calculated, fill = Target)) +
  ggplot2::stat_smooth(method = 'lm', lty = 11, lwd = .5, col = 'black', fill = 'grey50') +
  ggplot2::geom_point(shape = 21, size = 2) +
  ggplot2::scale_x_continuous(trans = scales::log10_trans(), limits = c(1, 10E8), expand = c(0,0)) +
  ggplot2::scale_y_continuous(trans = scales::log10_trans(), limits = c(1, 10E8), expand = c(0,0)) +
  ggplot2::scale_fill_manual(values = c('#FF1F5B', '#00CD6D', '#009ADE', '#FFC51E'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  ggplot2::labs(x = 'Predicted copy-number (log<sub>10</sub>)', y = 'Detected copy-number (log<sub>10</sub>)') +
  ggpmisc::stat_poly_eq(formula = x~y, rr.digits = 3, aes(label = ..rr.label..), parse = TRUE) +  
  ggplot2::facet_wrap(.~Target, nrow = 1) +
  ggplot2::theme(
    legend.position = 'bottom',
    legend.direction = 'horizontal',
    text = ggplot2::element_text(size = 9, family = 'Helvetica', face = 'bold'),
    axis.text.x = ggtext::element_markdown(),
    axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
    axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
    strip.text = ggtext::element_textbox_simple(width = NULL, halign = .5),
    panel.grid.major.x = ggplot2::element_line(colour = 'grey90', linetype = 'dotted'),
    panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
    panel.border = ggplot2::element_rect(fill = NA, colour = NA),
    strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.key = ggplot2::element_blank()
  )