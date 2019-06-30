# OM-GRC v2 ==============================================================================
# Code associated with the om-rgc v2 and tara prok metaT paper
# Figure X: Mapping rates ================================================================

rm(list = ls())
if (basename(getwd()) != 'analysis'){
  setwd('analysis')
}

# Libraries ------------------------------------------------------------------------------

library(tidyverse)
library(patchwork)
source('lib/Cell_Press_Guidelines.R')
source('reviewed_analysis/FigureX2A_Mapping_rates.R') # Load figure 2A
source('reviewed_analysis/FigureX2BC_Description_catalog.R') # Load figures 2B & 2C

# Build and save the figure --------------------------------------------------------------

figure_X2 = wrap_plots(figure_X2_A, figure_X2_B, figure_X2_C, widths = c(1, 2, 3))

ggsave(filename = '../results/figures/Figure_X2_Catalog_description.raw.pdf', plot_to_save, width = two_col, height = 90, units = col_unit)

