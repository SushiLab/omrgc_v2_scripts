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

# Global Variables -----------------------------------------------------------------------

mapping_rates = googlesheets::gs_read(ss = googlesheets::gs_title("aligment_stats"), skip = 4)
coding_perc = 87 # The average perc. of a genome that actually codes. Extracted from lit. and IMG.
# Matched metaG/metaT samples
matched_samples = readxl::read_xlsx("../results/paper_tables/submission/Table_S1.xlsx", sheet = "Table_W6")
plot_dest = '../results/figures/Figure_X2_A_mapping_rates.pdf'
table_dest = '../results/tables/Table_numbers_for_mapping_rates.tsv'
col_datasets = c("#D99233", "#4199DD")# MetaG, MetaT
vir_cols = viridis::viridis(2, begin = .2, end =.9)
vir_grey = c(DescTools::ColToGrey(vir_cols[length(vir_cols)]), vir_cols[1]) # Last color as greyscale

# Functions ------------------------------------------------------------------------------

check_table <- function(table){
  interest_sample_ids = Reduce(intersect, list(table$Sample,
                                               table$sample,
                                               table$Sample_1))
  if (length(interest_sample_ids) == nrow(table)) {
    print('All good to go')
  } else {
    warning('Tûût! You need to make sure your ids are comparable...
            for instance, did you remove META[G|T]/ from `Sample1`?')
  }
}


# Prepare mapping data -------------------------------------------------------------------

# Format Sample_1:
mapping_rates$Sample_1 = gsub('META[G|T]/', '', mapping_rates$Sample_1)
check_table(mapping_rates)

Formatted_table = 
  select(mapping_rates, c(Sample, `Fraction Aligned MARREF`)) %>%
  left_join(select(mapping_rates, c(sample, fraction_mapped_inserts, PANGAEA_ID)),
            by = c('Sample' = 'sample')) %>%
  left_join(select(mapping_rates, c(Sample_1, `Fraction Aligned DELMONT MAGS`)),
            by = c('Sample' = 'Sample_1')) %>%
  mutate(Dataset = ifelse(grepl('_G$', Sample), 'Metagenomes', 'Metatranscriptomes')) %>%
  # Factor in the coding fraction of the genomes in case of metagenomes to compare to catalog
  mutate(`MarRef\n(Ref. Genomes)` = ifelse(Dataset == 'Metagenomes',
                                           coding_perc/100*`Fraction Aligned MARREF`,
                                           `Fraction Aligned MARREF`),
         `Delmont et al. 2018\n(MAGs)` = ifelse(Dataset == 'Metagenomes',
                                                coding_perc/100*`Fraction Aligned DELMONT MAGS`,
                                                `Fraction Aligned DELMONT MAGS`),
         `This study\n(Gene Catalog)` = fraction_mapped_inserts*100) %>%
  select(-c(`Fraction Aligned MARREF`, fraction_mapped_inserts, `Fraction Aligned DELMONT MAGS`)) %>%
  filter(PANGAEA_ID %in% unlist(strsplit(matched_samples$Barcode, '-')))

assertthat::assert_that(length(unlist(strsplit(matched_samples$Barcode, '-'))) == 2*nrow(matched_samples),
                        msg = "Number of Pangea Id not equal to twice the number of pairs...")
assertthat::assert_that(nrow(Formatted_table) == 2*nrow(matched_samples),
                        msg = "Number of samples not equal to twice the number of pairs...")

Plot_table = Formatted_table %>%
  gather(key = 'data', value = 'aligned_perc' , -c(Dataset, Sample, PANGAEA_ID)) %>%
  mutate(data = factor(data, levels = c('MarRef\n(Ref. Genomes)',
                                        'Delmont et al. 2018\n(MAGs)',
                                        'This study\n(Gene Catalog)')))

# plot -----------------------------------------------------------------------------------

plot_mapping = ggplot(Plot_table) +
  geom_boxplot(aes(x = data, y = aligned_perc, fill = Dataset, colour = Dataset),
               alpha = .7, size = 0.75*size_converter, outlier.size = 2*size_converter) +
  theme_minimal() +
  theme_cell +
  ylim(0, 100) +
  ylab('Average percentage of reads aligned (%)') +
  scale_fill_manual(values = col_datasets) +
  scale_colour_manual(values = col_datasets) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(1, 1), 
        legend.justification = c(1, 1),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        axis.ticks = element_blank())

plot_mapping

Formatted_table %>%
  dplyr::rename(`MarRef (Ref. Genomes)` = `MarRef\n(Ref. Genomes)`,
                `This study (Gene Catalog)` = `This study\n(Gene Catalog)`,
                `Delmont et al. 2018 (MAGs)` = `Delmont et al. 2018\n(MAGs)`) %>%
  gather(key = Study, value = rate, -c(Sample, PANGAEA_ID, Dataset)) %>%
  group_by(Study, Dataset) %>%
  summarize(`Median mapping rate` = median(rate),
            `Mean mapping rate` = mean(rate)) %>%
  write_tsv(table_dest)

#ggsave(filename = plot_dest, plot_mapping)

source('paper_analysis/FigureX2BCD_Description_catalog.R')

plot_to_save = wrap_plots(plot_mapping, figure_X2_B, figure_X2_C, widths = c(1, 2, 3))

ggsave(filename = '../results/figures/Figure_X2_Catalog_description.raw.pdf', plot_to_save, width = two_col, height = 90, units = col_unit)

