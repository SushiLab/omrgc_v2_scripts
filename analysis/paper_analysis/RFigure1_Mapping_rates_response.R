# OM-GRC v2 ==============================================================================
# Code associated with the om-rgc v2 and tara prok metaT paper
# Figure X: Mapping rates ================================================================

rm(list = ls())
if (basename(getwd()) != 'analysis'){
  setwd('analysis')
}

# Libraries ------------------------------------------------------------------------------

library(tidyverse)

# Global Variables -----------------------------------------------------------------------

mapping_rates = googlesheets::gs_read(ss = googlesheets::gs_title("aligment_stats"), skip = 4)
almeida_mapping = round(73)
pasolli_mapping = round(87.51)
coding_perc = 87
plot_dest = '../results/figures/Figure_R1_mapping_rates.pdf'
vir_cols = viridis::viridis(2, begin = .2, end =.9)
vir_grey = c(DescTools::ColToGrey(vir_cols[length(vir_cols)]), vir_cols[1]) # Las color as greyscale

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
  mutate(Dataset = ifelse(grepl('_G$', Sample), 'MetaG', 'MetaT')) %>%
  group_by(Dataset) %>%
  summarise(`MarRef\n(Refs)` = mean(`Fraction Aligned MARREF`),
            `This Study\n(Gene Catalog)` = mean(fraction_mapped_inserts)*100,
            `Delmont et al. 2018\n(MAGs)` = mean(`Fraction Aligned DELMONT MAGS`)) %>%
  mutate(`Almeida et al. 2019\n(MAGs)` = c(almeida_mapping, NA),
         `Pasolli et al. 2019\n(Refs & MAGs)` = c(pasolli_mapping, NA)) %>%
  filter(Dataset == 'MetaG')

Plot_table = mutate(Formatted_table,
                    `This Study\n(Gene Catalog)` = `This Study\n(Gene Catalog)` +
                      `This Study\n(Gene Catalog)`/87*(100 - coding_perc),
                    Type = 'Non Coding') %>%
  bind_rows(mutate(Formatted_table,
                   `MarRef\n(Refs)` = coding_perc/100*`MarRef\n(Refs)`,
                   `Delmont et al. 2018\n(MAGs)` = coding_perc/100*`Delmont et al. 2018\n(MAGs)`,
                   `Almeida et al. 2019\n(MAGs)` = coding_perc/100*`Almeida et al. 2019\n(MAGs)`,
                   `Pasolli et al. 2019\n(Refs & MAGs)` = coding_perc/100*`Pasolli et al. 2019\n(Refs & MAGs)`,
                   Type = 'Coding Sequences')) %>%
  gather(key = 'data', value = 'aligned_perc' , -c(Dataset, Type)) %>%
  mutate(data = factor(data, levels = c('Pasolli et al. 2019\n(Refs & MAGs)',
                                        'Almeida et al. 2019\n(MAGs)',
                                        'MarRef\n(Refs)',
                                        'Delmont et al. 2018\n(MAGs)',
                                        'This Study\n(Gene Catalog)')))
Plot_table$aligned_perc[Plot_table$Type == 'Non Coding'] = 
  Plot_table$aligned_perc[Plot_table$Type == 'Non Coding'] -
  Plot_table$aligned_perc[Plot_table$Type == 'Coding Sequences']

Plot_table$Type = factor(Plot_table$Type, levels = c('Non Coding', 'Coding Sequences'))

# plot -----------------------------------------------------------------------------------

plot_mapping = ggplot(Plot_table) +
  geom_bar(aes(x = data, y = aligned_perc, fill = Type),
           stat = 'identity') +
  theme_bw() +
  ylim(0, 100) +
  annotate("segment", x = c(.75, 2.75), xend = c(2.25, 5.25), y = c(95, 95), yend = c(95, 95)) +
  annotate("text", x = c(1.5, 4), y = c(97.5, 97.5), label = c('Human Gut', 'Ocean')) +
  ylab('Average percentage of reads aligned (%)') +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) + 
  scale_fill_manual(values = vir_grey)
plot_mapping

ggsave(filename = plot_dest, plot_mapping, width = 8, height = 6, units = 'in')
