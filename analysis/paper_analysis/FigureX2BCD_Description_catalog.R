# OM-GRC v2 ==============================================================================
# Code associated with the om-rgc v2 and tara prok metaT paper
# Figure X2: Catalog descr. ==============================================================

#rm(list = ls()) to be sourced
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
                 "GC (Gene Clusters, 21.8 %)", "OG (eggNOG\nOrthologous Groups, 60.94 %)",
                 "KO (KEGG Orthologs, 23.59 %)")
pal<-sushi.palette()[c(14,1:13,15)]
vir_init = viridis(n = 10, direction = -1)
vir_init = c(DescTools::ColToGrey(vir_init[1]), vir_init[-1])
vir_6 = vir_init[c(1:5, 8)]
vir_4 = vir_init[c(6,7,9,10)]
mag_init = viridis(n = 4, direction = -1, option = "magma",
                   begin = .2, end = .8)
mag_5 = c(rep(DescTools::ColToGrey(vir_init[1]), 2), mag_init[-1])
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
  rbind(c("LUCA", round(prop_luca*sum(tax.stats.dom.raw$n)))) %>% # Add LUCA
  mutate(n = as.numeric(n),
         type = "Domain",
         Domain = ifelse(Domain == "", "No annotation", Domain)) %>%
  mutate(n = ifelse(Domain == "No annotation", n - round(prop_luca*sum(tax.stats.dom.raw$n)), n)) %>% # Remove LUCA annotation from Unannotated
  mutate(perc = 100*n/sum(n)) %>% 
  mutate(Domain = paste0(Domain, " (", round(perc,2), " %)"))
assertthat::assert_that(all(relevel_domain %in% tax.stats.dom$Domain))
tax.stats.dom = as_tibble(mutate(tax.stats.dom, Domain = fct_relevel(Domain, relevel_domain)))

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
  mutate(polar = ifelse(abs(Latitude) >= polarity, "polar", "non polar"))

# Generate plot --------------------------------------------------------------------------

p1<-ggplot(dades,aes(x=sample.num,y=n.genes,col=polar)) +
  geom_point(size = 1) +
  scale_color_manual(values=c(non_polar_col,polar_col)) +
  xlab("Prok. enriched samples") +
  ylab("Number of genes") +
  ylim(0,50000000) +
  geom_vline(xintercept = 139.5, linetype = 2, size = .75*size_converter) +
  annotate("segment", x = 40, xend = 180,y = 47000000, yend = 47000000,
           arrow = arrow(type = 'closed', angle = 20, length = unit(7, 'pt')), size = .75*size_converter) +
  annotate("text", x = 20, y = 47000000, label = "OM-RGC.v2\n(47M genes)", fontface = 'bold', size = 7*size_converter) +
  theme_bw() +
  theme_cell +
  theme(legend.position = c(0.3,0.2),
        legend.title = element_blank(),
        axis.text.y = element_text(angle = 90))

tibble_plot23 = rbind(dplyr::rename(tax.stats.dom, Taxon = Domain),
                      dplyr::rename(tax.stats.phyl, Taxon = Phylum)) %>%
  mutate(Taxon = factor(Taxon,
                        levels = c("Domain", levels(tax.stats.dom$Domain), " ",
                                   "Phylum (>2%)", levels(tax.stats.phyl$Phylum)))) # fix for the legend...

p23_main <- ggplot(tibble_plot23, aes(x="Dummy",y=perc,fill=Taxon)) +
  geom_bar(stat = "identity",position = position_stack()) +
  theme_minimal() +
  theme_cell +
  facet_wrap(~type) + # removed ", strip.position = 'bottom'" as it doesn't align well in combined plots
  scale_fill_manual(values=c(vir_6, vir_4)) +
  ylab("Percentage of genes") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = 'None')

p23_first_legend_pre <- tibble_plot23 %>%
  filter(Taxon %in% tax.stats.dom$Domain) %>%
  ggplot(aes("Dummy", fill = Taxon)) + geom_bar() + 
  scale_fill_manual(values = vir_6, name = "Domain") + theme_cell
p23_first_legend = ggdraw() + draw_plot(get_legend(p23_first_legend_pre))

p23_second_legend_pre <- tibble_plot23 %>%
  filter(Taxon %in% tax.stats.phyl$Phylum) %>%
  ggplot(aes("Dummy", fill = Taxon)) + geom_bar() + 
  scale_fill_manual(values = vir_4, name = "Phylum (>2%)") + theme_cell
p23_second_legend = ggdraw() + draw_plot(get_legend(p23_second_legend_pre))


tibble_plot4 = func.stats.merged %>%
  mutate(annotation = gsub('KO \\(', 'KO \\(KEGG Orthologs, ', annotation)) %>%
  mutate(annotation = gsub('NOG \\(', 'OG \\(eggNOG\nOrthologous Groups, ', annotation)) %>%
  mutate(annotation = gsub('GF \\(', 'GC \\(Gene Clusters, ', annotation)) %>%
  mutate(type = factor(type, levels = rev(unique(type))),
         annotation = fct_relevel(annotation, relevel_func))

p4_main <- ggplot(tibble_plot4, aes(x="Dummy",y=perc,fill=annotation)) +
  geom_bar(stat = "identity",position = position_stack()) +
  theme_minimal() +
  theme_cell + 
  facet_wrap(~type) + # removed ", strip.position = 'bottom'" as it doesn't align well in combined plots
  scale_fill_manual(values=mag_5) +
  ylab("Percentage of genes") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = 'None')

p4_first_legend_pre <- tibble_plot4 %>%
  filter(type == "OG + GC") %>%
  ggplot(aes(type, fill = annotation)) + geom_bar() + 
  scale_fill_manual(values = mag_5[2:4], name = 'OG + GC') + theme_cell
p4_first_legend = ggdraw() + draw_plot(get_legend(p4_first_legend_pre)) # using cowplot

p4_second_legend_pre <- tibble_plot4 %>%
  filter(type == "KO") %>%
  ggplot(aes(type, fill = annotation)) + geom_bar() + 
  scale_fill_manual(values =  mag_5[c(1,5)], name = 'KO') + theme_cell
p4_second_legend = ggdraw() + draw_plot(get_legend(p4_second_legend_pre))


figure_X2_B = p1 

figure_X2_C = 
  (p23_main | p4_main | (p23_first_legend / p23_second_legend / p4_first_legend / p4_second_legend)) +
  plot_layout(widths = c(1, 1, 2))

ggsave('../results/figures/Figure_X2_B_Catalog_accum.pdf', figure_X2_B)
ggsave('../results/figures/Figure_X2_C_Taxo_annot.pdf', figure_X2_C)
