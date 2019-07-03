# OM-GRC v2 ==============================================================================
# Code associated with the om-rgc v2 and tara prok metaT paper
# Figure X2: Catalog descr. ==============================================================

if (basename(getwd()) != 'analysis'){
  setwd('analysis')
}

# Libraries ------------------------------------------------------------------------------

library(data.table)
library(EcolUtils)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(cowplot)
library(viridis)
source("lib/varpart.sqr.euc_functions.R")
source("lib/Cell_Press_Guidelines.R")
source("lib/sushipal.R")

# Variables ------------------------------------------------------------------------------

polarity = 53 # Above 53 degrees north, stations are considered polar: starts with st155
prop_luca = 0.0336 # Proportion of tax annotation to LUCA (3.36 %). Extracted from the figure
# as it's not in the tables... (FIXME?)
# Order for plot:
relevel_domain = c("No annotation (27.48 %)",
                   "Archaea (2.03 %)",
                   "Eukaryota (2.61 %)",
                   "LUCA (3.36 %)",
                   "Viruses (5.26 %)",
                   "Bacteria (59.26 %)")
relevel_phylum = "No annotation (52.06 %)" # extracted from the plot
relevel_func = c("No annotation (76.41 %)","No annotation (17.26 %)",
                 "GC (Gene Clusters, 21.8 %)", "OG (eggNOG 60.94 %)",
                 "KO (KEGG 23.59 %)")
pal<-sushi.palette()
ggplot() + geom_point(aes(x=1:15,y=1:15,col=as.factor(1:15)),size=5)+
  scale_color_manual(values=pal)+theme_bw()+theme(legend.position = 'none')
vir_init = viridis(n = 10, direction = -1)
init_grey = DescTools::ColToGrey(vir_init[1])
vir_init = c(init_grey, vir_init[-1])
vir_6 = vir_init[c(1:5, 8)]
pal_6 = pal[c(14, 4, 5, 3, 2, 1)]
# colors from 2015 paper c('#D9D9D9', '#CDE5C4', '#F9CDE0', '#BEBBDB', '#FAF7B5', '#4CB049')
# get darker for the pale ones
col_6 = c('#D9D9D9', '#c7f2b7', '#F9CDE0', '#BEBBDB', '#FAF7B5', '#4CB049')
vir_4 = vir_init[c(6,7,9,10)]
pal_4 = pal[c(1, 2, 3, 4)]
mag_init = viridis(n = 4, direction = -1, option = "magma", begin = .2, end = .8)
mag_5 = c(rep(init_grey, 2), mag_init[-1])
pal_5 = pal[c(14, 14, 1, 3, 2)]
col_5 = c('#D9D9D9', '#D9D9D9', pal[9], pal[2], pal[1])
polar_col = "#136FBA"
non_polar_col = "#F98B04"


# Load data ------------------------------------------------------------------------------

# functional annotations
func.stats<-fread("../data/processed/stats.func.annotation.tsv",header=T,sep="\t",data.table = F)
# taxo annotations
tax.stats.dom.raw<-fread("../data/processed/stats.tax.annotation.Domain.tsv",header=T,sep="\t",data.table = F)
tax.stats.phyl.raw<-fread("../data/processed/stats.tax.annotation.Phylum.tsv",header=T,sep="\t",data.table = F)
tax.stats.class.raw<-fread("../data/processed/stats.tax.annotation.Class.tsv",header=T,sep="\t",data.table = F)
# gene accumulation
dades.raw<-fread("zcat < ../data/processed/stats.accum.genes.tsv.gz",sep="\t",header=T,data.table = F)
# metadata
env.mat<-fread("zcat < ../data/processed/KO_env.mat.metaG.txt.gz",sep="\t", header=T)

# Prepare data ---------------------------------------------------------------------------

# Functional annotations ####
func.stats.ko<-func.stats[c(1,2),]
func.stats.ko[1,2]<-func.stats.ko[1,2]-func.stats.ko[2,2]
func.stats.ko[1,1]<-"No annotation"
func.stats.ko$type<-"KO"
func.stats.ko<-func.stats.ko %>% mutate(perc=100*n/sum(n))

func.stats.egg<-func.stats[c(1,3,4),]
func.stats.egg[1,2]<-func.stats.egg[1,2]-sum(func.stats.egg[c(2,3),2])
func.stats.egg[1,1]<-"No annotation"
func.stats.egg$type<-"OG + GC"
func.stats.egg<-func.stats.egg %>% mutate(perc=100*n/sum(n))

func.stats.merged<-rbind(func.stats.ko,func.stats.egg)
func.stats.merged<-func.stats.merged %>% mutate(annotation=paste(annotation," (",round(perc,2)," %)",sep=""))

# Taxo Domain level ####
tax.stats.dom = tax.stats.dom.raw %>%
  rename(Taxonomy = Domain) %>%
  rbind(c("LUCA", round(prop_luca*sum(tax.stats.dom.raw$n)))) %>% # Add LUCA
  mutate(n = as.numeric(n),
         type = "Taxonomy",
         Taxonomy = ifelse(Taxonomy == "", "No annotation", Taxonomy)) %>%
  mutate(n = ifelse(Taxonomy == "No annotation", n - round(prop_luca*sum(tax.stats.dom.raw$n)), n)) %>% # Remove LUCA annotation from Unannotated
  mutate(perc = 100*n/sum(n)) %>% 
  mutate(Taxonomy = paste0(Taxonomy, " (", round(perc,2), " %)"))
assertthat::assert_that(all(relevel_domain %in% tax.stats.dom$Taxonomy))
tax.stats.dom = as_tibble(mutate(tax.stats.dom, Taxonomy = fct_relevel(Taxonomy, relevel_domain)))

# Taxo Phylum/Class level ####
tax.stats.class = tax.stats.class.raw %>%
  dplyr::rename(Phylum = Class) # rename class for merge with phylum
colnames(tax.stats.class)[1]<-"Phylum"

tax.stats.phyl = as_tibble(tax.stats.phyl.raw) %>%
  filter(Phylum != "Proteobacteria") %>%
  rbind(filter(tax.stats.class, grepl("proteobacteria", Phylum))) %>%
  mutate(n = ifelse(Phylum == "", n + # Add the Proteobacteria not classified at the class level as unknown
                      sum(filter(tax.stats.phyl.raw, Phylum == 'Proteobacteria')$n) - sum(filter(tax.stats.class, grepl('proteobacteria', Phylum))$n), n),
         type = "Phylum") %>%
  mutate(Phylum = ifelse(Phylum == "", "No annotation", Phylum)) %>%
  #mutate(Phylum = ifelse(n < 200000, "Other", Phylum)) %>% # Initial plot
  group_by(Phylum, type) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  mutate(perc = 100*n/sum(n)) %>%
  mutate(Phylum = paste0(Phylum," (",round(perc,2)," %)")) %>%
  filter(perc >= 2) %>% # For new plot
  mutate(Phylum = fct_reorder(Phylum, dplyr::desc(n)))
assertthat::assert_that(all(relevel_phylum %in% tax.stats.phyl$Phylum)) # FIXME: Expects (52.06 %)
tax.stats.phyl = mutate(tax.stats.phyl, Phylum = fct_relevel(Phylum, relevel_phylum)) 
tax.stats.phyl = filter(tax.stats.phyl, Phylum != relevel_phylum) # For new plot

# Gene accumulation ####
dades = dades.raw %>%
  mutate(sample.num = 1:nrow(dades.raw)) %>%
  cbind(env.mat) %>%
  mutate(polar = ifelse(abs(Latitude) >= polarity, "Polar", "Non-polar"))

# Format tables for plots
tibble_plot23 = rbind(dplyr::rename(tax.stats.dom, Taxon = Taxonomy),
                      dplyr::rename(tax.stats.phyl, Taxon = Phylum))

tibble_plot4 = func.stats.merged %>%
  mutate(annotation = gsub('KO \\(', 'KO \\(KEGG ', annotation)) %>%
  mutate(annotation = gsub('NOG \\(', 'OG \\(eggNOG ', annotation)) %>%
  mutate(annotation = gsub('GF \\(', 'GC \\(Gene Clusters, ', annotation)) %>%
  mutate(type = factor(type, levels = rev(unique(type))),
         annotation = fct_relevel(annotation, relevel_func))

tibble_plot234 = rbind(rename(tibble_plot23, Legend = Taxon), rename(tibble_plot4, Legend = annotation)) %>%
  mutate(type = factor(type, levels = c("Taxonomy", "Phylum", "OG + GC", "KO"))) %>%
  filter(type != 'Phylum') # Remove phylum

# Generate plot --------------------------------------------------------------------------

p1<-ggplot(dades,aes(x=sample.num,y=n.genes,col=polar)) +
  geom_point(size = 2*size_converter) +
  scale_color_manual(values=c(non_polar_col,polar_col), guide = guide_legend(label.position = 'left')) +
  xlab("Prokaryote-enriched samples") +
  ylab(expression(paste('Number of genes (',10^{6},')'))) +
  geom_vline(xintercept = 139.5, linetype = 2, size = .5*size_converter) +
  annotate("segment", x = 70, xend = 180,y = 47000000, yend = 47000000,
           arrow = arrow(type = 'closed', angle = 20, length = unit(7, 'pt')), size = .75*size_converter) +
  annotate("text", x = 35, y = 47000000, label = "OM-RGC.v2\n(47M genes)", fontface = 'bold', size = 7*size_converter) +
  scale_y_continuous(limits = c(0,47000000), expand = expand_scale(mult = c(0, .05)),
                     breaks = c(0, 10000000, 20000000, 30000000, 40000000, 47000000),
                     labels = c('0', '10', '20', '30', '40', '47')) +
  theme_minimal() +
  theme_cell +
  theme(legend.position = c(0.3,0.2),
        legend.title = element_blank(),
        legend.background = element_rect(colour = 'black', fill = 'white', size = .5*size_converter),
        axis.title.y = element_text(size = unit(9, 'pt')),
        axis.line = element_line(colour = 'black', size = 0.5*size_converter))
p1

p234_main <- ggplot(tibble_plot234, aes(x=type,y=perc,fill=Legend)) +
  geom_bar(stat = "identity",position = position_stack()) +
  theme_minimal() +
  theme_cell +
  #scale_fill_manual(values=c(pal_6, pal_4, pal_5)) +
  scale_fill_manual(values=c(col_6, col_5)) +
  scale_y_continuous(limits = c(0, 100), expand = expand_scale(mult = c(0, .05))) +
  ylab("Percentage of genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = unit(9, 'pt')),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = unit(9, 'pt')),
        axis.line.x = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = 'None')
p234_main


p23_first_legend_pre <- tibble_plot23 %>%
  filter(Taxon %in% tax.stats.dom$Taxonomy) %>%
  ggplot(aes("Dummy", fill = Taxon)) + geom_bar() + 
  scale_fill_manual(values = col_6, name = "Taxonomy") + theme_cell
p23_first_legend = ggdraw() + draw_plot(get_legend(p23_first_legend_pre))

p23_second_legend_pre <- tibble_plot23 %>%
  filter(Taxon %in% tax.stats.phyl$Phylum) %>%
  ggplot(aes("Dummy", fill = Taxon)) + geom_bar() + 
  scale_fill_manual(values = pal_4, name = "Phylum (>2%)") + theme_cell
p23_second_legend = ggdraw() + draw_plot(get_legend(p23_second_legend_pre))

p4_first_legend_pre <- tibble_plot4 %>%
  filter(type == "OG + GC") %>%
  ggplot(aes(type, fill = annotation)) + geom_bar() + 
  scale_fill_manual(values = col_5[2:4], name = 'OG + GC') + theme_cell
p4_first_legend = ggdraw() + draw_plot(get_legend(p4_first_legend_pre)) # using cowplot

p4_second_legend_pre <- tibble_plot4 %>%
  filter(type == "KO") %>%
  ggplot(aes(type, fill = annotation)) + geom_bar() + 
  scale_fill_manual(values =  col_5[c(1,5)], name = 'KO') + theme_cell
p4_second_legend = ggdraw() + draw_plot(get_legend(p4_second_legend_pre))

# Get everything organised:
figure_X2_B = p1 
figure_X2_C = 
  # (p234_main | (p23_first_legend / p23_second_legend / p4_first_legend / p4_second_legend) + 
  #    plot_layout(height = c(6, 4, 4, 2))) +
  (p234_main | (p23_first_legend / p4_first_legend / p4_second_legend / plot_spacer()) + 
     plot_layout(height = c(6, 4, 3, 2))) +
  plot_layout(widths = c(2, 4))

