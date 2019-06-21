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

# Global Variables -----------------------------------------------------------------------

mapping_rates = googlesheets::gs_read(ss = googlesheets::gs_title("aligment_stats"), skip = 4)
coding_perc = 87
plot_dest = '../results/figures/Figure_X2_A_mapping_rates.pdf'
table_dest = '../results/tables/Table_numbers_for_mapping_rates.tsv'
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
  left_join(select(mapping_rates, c(sample, fraction_mapped_inserts)),
            by = c('Sample' = 'sample')) %>%
  left_join(select(mapping_rates, c(Sample_1, `Fraction Aligned DELMONT MAGS`)),
            by = c('Sample' = 'Sample_1')) %>%
  mutate(Dataset = ifelse(grepl('_G$', Sample), 'Metagenomes', 'Metatranscriptomes')) %>%
  group_by(Dataset) %>%
  summarise(`MarRef\n(Ref. Genomes)` = median(`Fraction Aligned MARREF`),
            `This study\n(Gene Catalog)` = median(fraction_mapped_inserts)*100,
            `Delmont et al. 2018\n(MAGs)` = median(`Fraction Aligned DELMONT MAGS`))
Formatted_table = Formatted_table %>%
  mutate(`This study\n(Gene Catalog)` = `This study\n(Gene Catalog)` +
           `This study\n(Gene Catalog)`/87*(100 - coding_perc),
         Type = 'Non Coding') %>%
  bind_rows(mutate(Formatted_table,
                   `MarRef\n(Ref. Genomes)` = coding_perc/100*`MarRef\n(Ref. Genomes)`,
                   `Delmont et al. 2018\n(MAGs)` = coding_perc/100*`Delmont et al. 2018\n(MAGs)`,
                   Type = 'Coding Sequences'))

Plot_table = Formatted_table %>%
  gather(key = 'data', value = 'aligned_perc' , -c(Dataset, Type)) %>%
  mutate(data = factor(data, levels = c('MarRef\n(Ref. Genomes)',
                                        'Delmont et al. 2018\n(MAGs)',
                                        'This study\n(Gene Catalog)')))

Plot_table$aligned_perc[Plot_table$Type == 'Non Coding'] = 
  Plot_table$aligned_perc[Plot_table$Type == 'Non Coding'] -
  Plot_table$aligned_perc[Plot_table$Type == 'Coding Sequences']

Plot_table$Type = factor(Plot_table$Type, levels = c('Non Coding', 'Coding Sequences'))

Plot_table = Plot_table %>%
  filter(!(Dataset == 'Metatranscriptomes' & Type == 'Non Coding'))

# plot -----------------------------------------------------------------------------------

plot_mapping = ggplot(Plot_table) +
  geom_bar(aes(x = data, y = aligned_perc, fill = Type),
           stat = 'identity') +
  theme_bw() +
  ylim(0, 100) +
  ylab('Average percentage of reads aligned (%)') +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) + 
  scale_fill_manual(values = vir_grey) +
  scale_colour_manual(values = c('black', NA)) +
  facet_wrap(~Dataset) +
  theme(legend.position = c(.95, .95), 
        legend.justification = c(1, 1),
        axis.text.x = element_text(angle = 35, hjust = 1))

plot_mapping

Formatted_table %>%
  dplyr::rename(`MarRef (Ref. Genomes)` = `MarRef\n(Ref. Genomes)`,
                `This study (Gene Catalog)` = `This study\n(Gene Catalog)`,
                `Delmont et al. 2018 (MAGs)` = `Delmont et al. 2018\n(MAGs)`) %>%
  write_tsv(table_dest)
ggsave(filename = plot_dest, plot_mapping)

source('paper_analysis/FigureX2BCD_Description_catalog.R')

plot_to_save = wrap_plots(plot_mapping, figure_X2_B) /
  wrap_plots(figure_X2_C, figure_X2_D)


ggsave(filename = '../results/figures/Figure_X2_Catalog_description.raw.pdf', plot_to_save, width = 12, height = 10, units = 'in')

