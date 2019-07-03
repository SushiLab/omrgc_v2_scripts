# OM-GRC v2 ==============================================================================
# Code associated with the om-rgc v2 and tara prok metaT paper
# Figure 4: Distance partitioning ========================================================

rm(list = ls())
if (basename(getwd()) != 'analysis'){
  setwd('analysis')
}

# Libraries ------------------------------------------------------------------------------

library(geosphere)
library(tidyverse)
library(data.table)
library(patchwork)
library(cowplot)
source("lib/sushipal.R")
source("lib/varpart.sqr.euc_functions.R")
source("lib/Cell_Press_Guidelines.R")
palette(sushi.palette(alpha=0.7)[c(2,3,4,1,14)])

# Variables ------------------------------------------------------------------------------

polar_col = "#136FBA"
non_polar_col = "#F98B04"
table_dest = '../results/tables/Table_numbers_for_bins_temp.tsv'

# Load data ------------------------------------------------------------------------------

metaG.norm.match.log2<-fread("zcat < ../data/processed/NOG_metaG.norm.match.log2.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match.log2)<-metaG.norm.match.log2$V1
metaG.norm.match.log2<-metaG.norm.match.log2[,-1]

metaG.norm.match<-fread("zcat < ../data/processed/NOG_metaG.norm.match.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaG.norm.match)<-metaG.norm.match$V1
metaG.norm.match<-metaG.norm.match[,-1]

metaT.norm.match.log2<-fread("zcat < ../data/processed/NOG_metaT.norm.match.log2.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match.log2)<-metaT.norm.match.log2$V1
metaT.norm.match.log2<-metaT.norm.match.log2[,-1]

metaT.norm.match<-fread("zcat < ../data/processed/NOG_metaT.norm.match.txt.gz",header=T,sep="\t",data.table = F)
rownames(metaT.norm.match)<-metaT.norm.match$V1
metaT.norm.match<-metaT.norm.match[,-1]

ratio.mat<-fread("zcat < ../data/processed/NOG_ratio.mat.txt.gz",header=T,sep="\t",data.table = F)
rownames(ratio.mat)<-ratio.mat$V1
ratio.mat<-ratio.mat[,-1]

env.mat.match<-fread("zcat < ../data/processed/NOG_env.mat.match.txt.gz",header=T,sep="\t",data.table = F,stringsAsFactors = T)
rownames(env.mat.match)<-env.mat.match$V1
env.mat.match<-env.mat.match[,-1]
env.mat.match$epi<-env.mat.match$Layer
levels(env.mat.match$epi)<-c("EPI","MES","EPI","EPI")
env.mat.match$date<-as.POSIXct(substr(as.character(env.mat.match$Event.date),1,10))
env.mat.match$doy<-yday(env.mat.match$date)
env.mat.match$daylength<-as.numeric(daylength(lat = env.mat.match$Latitude,doy=env.mat.match$doy))
env.mat.match$hour<-as.numeric(substr(env.mat.match$Event.date,12,13))
env.mat.match<-env.mat.match[match(rownames(metaG.norm.match),env.mat.match$Barcode),]

# EPI - Variance partitioning: Sample bins through temperature (each layers) -------------

# Compute the distances, takes time
res.epi<-varpart.sqr.euc.all(metaT.norm.match.log2[env.mat.match$Layer %in% c("SRF","DCM"),],
                             metaG.norm.match.log2[env.mat.match$Layer %in% c("SRF","DCM"),],
                             ratio.mat[env.mat.match$Layer %in% c("SRF","DCM"),])
# Format the results
res.epi.comp <- as_tibble(res.epi$components) %>%
  mutate(comparison_id = paste0(sample1, ':', sample2)) %>%
  gather("Component","value",-sample1,-sample2,-comparison_id) %>%
  mutate(comparison_id_w_component = paste0(comparison_id, ':', Component),
         Layer = 'EPI')

# Bin the samples by temperature
mat.env<-env.mat.match[env.mat.match$Layer %in% c("SRF","DCM"),]
n<-15
ordre<-order(mat.env$Temperature)
res.epi.binned <- NULL
for (i in 1:(nrow(mat.env)-n)){
  smpls<-mat.env$Barcode[ordre[i:(i+n-1)]]
  pos<-which(res.epi.comp$sample1 %in% smpls & res.epi.comp$sample2 %in% smpls)
  tmp<-res.epi.comp[pos,]
  tmp$bin<-i
  tmp$median.temp<-median(mat.env$Temperature[ordre[i:(i+n-1)]])
  tmp$width<-max(mat.env$Temperature[ordre[i:(i+n-1)]])-min(mat.env$Temperature[ordre[i:(i+n-1)]])
  res.epi.binned<-rbind(res.epi.binned,tmp)
}

# Additional formatting for downstream analyses
add_polarity <- function(sample1, sample2, polar1, polar2){
  n = length(unique(c(sample1, sample2)))
  assertthat::assert_that(n == 15,
                          msg = 'Bins should be of size 15...')
  select_samples = !duplicated(c(sample1, sample2))
  init_polarity = c(polar1, polar2)[select_samples]
  polarity = sum(init_polarity == 'polar') / 15
  return(polarity)
}

res.epi.binned = res.epi.binned %>%
  mutate(polar1 = as.character(env.mat.match$polar[match(sample1,env.mat.match$Barcode)]),
         polar2 = as.character(env.mat.match$polar[match(sample2,env.mat.match$Barcode)])) %>%
  mutate(polar = paste0(polar1, '.', polar2)) %>%
  group_by(bin) %>%
  mutate(polarity = add_polarity(sample1, sample2, polar1, polar2)) %>%
  ungroup()

assertthat::validate_that(nrow(res.epi.binned) == length(unique(res.epi.binned$comparison_id_w_component)),
                          msg = 'Note that some samples are compared to each other several times.')

# Look at an example:
# res.epi.binned %>%
#   filter(comparison_id == "TARA_B110000971-TARA_B110000969:TARA_B100000787-TARA_B100000787") %>%
#   View()

# Plot distance distributions (violin plot) ----------------------------------------------

tibble_violin_plot = res.epi.binned %>% 
  filter(!duplicated(comparison_id_w_component)) %>%
  filter(polar %in% c("polar.polar","non polar.non polar")) %>%
  filter(Component != 'interaction') %>%
  mutate(polar = factor(polar, levels = c("polar.polar","non polar.non polar"),
                        labels = c('Polar', 'Non-polar')),
         Component = factor(Component, labels = c('Turnover', 'Acclimatization')))

assertthat::assert_that(nrow(tibble_violin_plot) < 2*length(unique(res.epi.binned$comparison_id)),
                        msg = 'Should have less than twice the number of comparison as you remove interaction and non homogeneous comparisons.')

t.test(filter(tibble_violin_plot, Component == 'Turnover' & polar == 'Polar')$value,
       filter(tibble_violin_plot, Component == 'Turnover' & polar == 'Non-polar')$value)
wilcox.test(filter(tibble_violin_plot, Component == 'Turnover' & polar == 'Polar')$value,
            filter(tibble_violin_plot, Component == 'Turnover' & polar == 'Non-polar')$value)
t.test(filter(tibble_violin_plot, Component == 'Acclimatization' & polar == 'Polar')$value,
       filter(tibble_violin_plot, Component == 'Acclimatization' & polar == 'Non-polar')$value)
wilcox.test(filter(tibble_violin_plot, Component == 'Acclimatization' & polar == 'Polar')$value,
            filter(tibble_violin_plot, Component == 'Acclimatization' & polar == 'Non-polar')$value)

# Also check for interaction
t.test(filter(res.epi.binned, !duplicated(comparison_id_w_component) & Component == 'interaction' & polar == 'polar.polar')$value,
       filter(res.epi.binned, !duplicated(comparison_id_w_component) & Component == 'interaction' & polar == 'non polar.non polar')$value)

p_violin = ggplot(tibble_violin_plot) +
  geom_violin(aes(x=polar,y=value,fill=polar),
              scale='width',draw_quantiles = c(0.5), size = .5*size_converter) +
  annotate("segment", x = 1, xend = 2, y = 6, yend = 6, size = .5*size_converter) +
  annotate("text", label = '***', x = 1.5, y = 6.05, size = 7*size_converter) +
  ylim(0,6.1) +
  scale_fill_manual(values = c(polar_col, non_polar_col)) +
  facet_wrap(~Component) +
  ylab('Distances') +
  theme_bw() +
  theme_cell +
  theme(plot.background = element_rect(colour = 'black', fill = 'white', size = .5*size_converter),
        panel.background = element_rect(size = .5*size_converter),
        panel.border = element_rect(size = .5*size_converter),
        strip.background = element_rect(size = .5*size_converter),
        strip.text.x = element_text(size = unit(6, 'pt'), margin = margin(2, 2, 2, 2, unit = 'pt')),
        legend.position = 'None',
        axis.title.y = element_text(size = unit(6, 'pt')),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1.3, size = unit(6, 'pt')),
        axis.ticks.x = element_blank(),
        plot.margin = margin(3, 3, 3, 3, unit = 'pt'))

p_violin

# Prepare the data for bin plot ----------------------------------------------------------

res.epi.binned.med <-res.epi.binned %>%
  group_by(bin,Component) %>%
  summarise(mean=mean(value),sd=sd(value),median=median(value),q1=quantile(value,probs = 0.25),q3=quantile(value,probs = 0.75)) %>%
  left_join(res.epi.binned)

wilcox_wrapper <- function(x, table = res.epi.binned.med){
  test <- wilcox.test(table$value[table$bin==x & table$Component=="Abundance"],
                      table$value[table$bin==x & table$Component=="Expression"])
  test$p.value
}

# Per bin p-values
pval.df <- tibble(bin = unique(res.epi.binned.med$bin),
                  p.val = sapply(unique(res.epi.binned.med$bin), wilcox_wrapper))

# FIXME: Wait.. why bonferroni + <0.01 ? Sounds very conservative. Replaced by Holm/0.05
pval.df = pval.df %>%
  mutate(p.val.holm =  p.adjust(p.val, method = "holm"),
         p.val.bonf =  p.adjust(p.val, method = "bonferroni")) %>%
  mutate(sign = ifelse(p.val.holm <= 0.05, "Significant\n(p < 0.05)", "Not Significant"))

res.epi.binned.med.simpl <- left_join(res.epi.binned.med, pval.df, by="bin") %>%
  dplyr::select(bin, Component, mean, sd, median, q1, q3, Layer, median.temp, width, p.val, p.val.holm, sign, polarity) %>%
  unique()

# FIXME: Why would interactions always be significant ??
res.epi.binned.med.simpl = mutate(res.epi.binned.med.simpl, sign = ifelse(Component == 'interaction', 'Significant', sign))

# Plot distance partitioning without interactions ----------------------------------------

tibble_ratio_plot = res.epi.binned.med.simpl %>% 
  filter(Component != "interaction") %>%
  group_by(bin) %>%
  summarize(Bin_temperature = unique(median.temp),
            Diff = 
              unique(median[Component == 'Abundance']) -
              unique(median[Component == 'Expression']),
            Ratio = 
              unique(median[Component == 'Abundance']) /
              unique(median[Component == 'Expression']),
            Significance = unique(sign),
            Polarity = unique(polarity))

# FIXME add color of the line/points based on % of polar samples.
p_ratio<-ggplot(tibble_ratio_plot) +
  geom_hline(yintercept = 1, linetype = 2, color = 'grey75', size = .75*size_converter) +
  geom_point(aes(x = Bin_temperature, y = Ratio, shape = Significance, fill = Polarity),
             size = 7*size_converter, stroke = .4*size_converter) +
  theme_minimal() +
  ylab("Abundance-based distance / Expression-based distance") +
  xlab("Median temperature of the bin (°C)") +
  scale_shape_manual(values = c(21,23)) +
  scale_fill_gradient2(low = non_polar_col, high = polar_col, mid = 'lightgrey', midpoint = 0.5,
                       breaks = c(0, 1), labels = c('All non-polar', 'All polar')) +
  labs(fill = 'Proportion of samples', shape = element_blank()) +
  theme_cell +
  theme(legend.position = c(.95, 1), 
        legend.justification = c(1, 1), 
        legend.key.height = unit(4, 'mm'),
        legend.key.width = unit(6, 'mm'),
        legend.box.background = element_rect(colour = "black", fill = "white", size = .4*size_converter),
        plot.margin = unit(c(1,2,1,1), 'lines'),
        axis.title = element_text(size = unit(7, 'pt')),
        axis.line = element_line(colour = 'black', size = 0.5*size_converter)) +
  annotate("segment", x = 30, xend = 30, y = 1.05, yend = 1.4, arrow = arrow(type = 'closed', angle = 20, length = unit(7, 'pt'))) +
  annotate("text", x = 31.75, y = 1.225, label = 'Turnover', angle = -90, size = 8*size_converter, ) +
  annotate("text", x = 30.75, y = 1.225, label = 'dominates', angle = -90, size = 8*size_converter, ) +
  annotate("segment", x = 30, xend = 30, y = .95, yend = .6, arrow = arrow(type = 'closed', angle = 20, length = unit(7, 'pt'))) +
  annotate("text", x = 31.75, y = .775, label = 'Acclimatization', angle = -90, size = 8*size_converter) +
  annotate("text", x = 30.75, y = .775, label = 'dominates', angle = -90, size = 8*size_converter) +
  coord_cartesian(xlim = c(0, 28), clip = 'off')

p_ratio

plot.with.inset <-
  ggdraw() +
  draw_plot(p_ratio) +
  draw_plot(p_violin, x = .13, y = .13, width = .35, height = .45)

plot.with.inset

ggsave("../results/figures/Figure_X6_distance_partitioning.raw.pdf", plot.with.inset, height = 90, width = one_half_col, unit = col_unit)


# Table to save results in main text -----------------------------------------------------

res.epi.binned.med.simpl %>% 
  ungroup() %>%
  filter(Component != "interaction") %>%
  summarize(Description = 'The width of a bin represents the max temp. diff. between two samples',
            `Min temperature diff` = min(width),
            `Median temperature diff` = median(width),
            `Max temperature diff` = max(width)) %>%
  write_tsv(table_dest)
