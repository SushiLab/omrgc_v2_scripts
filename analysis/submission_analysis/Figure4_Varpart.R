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
source("lib/sushipal.R")
source("lib/varpart.sqr.euc_functions.R")
palette(sushi.palette(alpha=0.7)[c(2,3,4,1,14)])

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

# EPI
res.epi<-varpart.sqr.euc.all(metaT.norm.match.log2[env.mat.match$Layer %in% c("SRF","DCM"),],
                             metaG.norm.match.log2[env.mat.match$Layer %in% c("SRF","DCM"),],
                             ratio.mat[env.mat.match$Layer %in% c("SRF","DCM"),])
res.epi.comp.norm<-res.epi$components %>% gather("Component","value",-sample1,-sample2)
res.epi.comp.norm$Layer<-"EPI"

mat.env<-env.mat.match[env.mat.match$Layer %in% c("SRF","DCM"),]
n<-15
ordre<-order(mat.env$Temperature)
res<-NULL
for (i in 1:(nrow(mat.env)-n)){
  smpls<-mat.env$Barcode[ordre[i:(i+n-1)]]
  pos<-which(res.epi.comp.norm$sample1 %in% smpls & res.epi.comp.norm$sample2 %in% smpls)
  tmp<-res.epi.comp.norm[pos,]
  tmp$bin<-i
  tmp$median.temp<-median(mat.env$Temperature[ordre[i:(i+n-1)]])
  tmp$width<-max(mat.env$Temperature[ordre[i:(i+n-1)]])-min(mat.env$Temperature[ordre[i:(i+n-1)]])
  res<-rbind(res,tmp)
}

# Plot distance distributions ------------------------------------------------------------

res$polar1<-env.mat.match$polar[match(res$sample1,env.mat.match$Barcode)]
res$polar2<-env.mat.match$polar[match(res$sample2,env.mat.match$Barcode)]
res$polar<-interaction(res$polar1,res$polar2)

ggplot(data=res %>% filter(polar %in% c("polar.polar","non polar.non polar")),aes(x=Component,y=value,fill=polar)) +
  geom_violin(scale = "width",draw_quantiles = c(0.025,0.5,0.975)) +
  theme_bw()

t.test(res$value[res$Component=="Abundance" & res$polar=="polar.polar"],res$value[res$Component=="Abundance" & res$polar=="non polar.non polar"])
t.test(res$value[res$Component=="Expression" & res$polar=="polar.polar"],res$value[res$Component=="Expression" & res$polar=="non polar.non polar"])
t.test(res$value[res$Component=="interaction" & res$polar=="polar.polar"],res$value[res$Component=="interaction" & res$polar=="non polar.non polar"])

# Plot distance partition with interactions ----------------------------------------------

res.med<-res %>%
  group_by(bin,Component) %>%
  summarise(mean=mean(value),sd=sd(value),median=median(value),q1=quantile(value,probs = 0.25),q3=quantile(value,probs = 0.75)) %>%
  as.data.frame() %>%
  left_join(res,by="bin") %>%
  dplyr::select(-Component.y)

fun<-function(x){
  test<-wilcox.test(res$value[res$bin==x & res$Component=="Abundance"],res$value[res$bin==x & res$Component=="Expression"])
  test$p.value
}
pval.df<-data.frame(bin=unique(res.med$bin),p.val=sapply(unique(res.med$bin),fun))
pval.df$p.val.fdr<-p.adjust(pval.df$p.val,method = "bonferroni")
pval.df$sign<-sapply(pval.df$p.val.fdr,function(x){if (x<0.01) "P<0.05" else "NS"})

res.med<-left_join(res.med,pval.df,by="bin") %>%
  dplyr::select(-sample1,-sample2,-value) %>%
  unique()
res.med$sign[res.med$Component.x=="interaction"]<-"P<0.05"

g1<-ggplot() +
  geom_point(data=res.med,aes(x=median.temp,y=median,col=Component.x,shape=sign)) +
  geom_smooth(data=res.med,aes(x=median.temp,y=median,col=Component.x),method="loess",span=0.1,se=F) +
  geom_ribbon(data=res.med,aes(x=median.temp,ymin=q1,ymax=q3,fill=Component.x),alpha=0.3) +
  #geom_text(data=res.med,aes(x=median.temp,y=3.2,label=sign),size=4) +
  theme_bw() +
  ylab("Component value") +
  xlab("Bins' median temperature") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c("#eaad53","#90e25a","gray")) +
  scale_fill_manual(values = c("#eaad53","#90e25a","gray")) +
  scale_shape_manual(values = c(3,19))

g1.later<-g1 +
  scale_y_continuous(limits = c(-4,4)) +
  labs(title = "Epipelagic")


# Plots for panels A and B: Oceans, Latitude, Temperature... -----------------------------

res.env<-res %>%
  dplyr::select(sample1,sample2,bin,median.temp) %>%
  gather(key = "a",value="Barcode",-bin,-median.temp) %>%
  dplyr::select(-a) %>%
  left_join(env.mat.match,by="Barcode") %>%
  unique()

res.lat<-res.env %>%
  group_by(bin,median.temp) %>%
  summarise(median.lat=median(abs(Latitude)),q1=quantile(abs(Latitude),probs = 0.25),q3=quantile(abs(Latitude),probs = 0.75)) %>%
  as.data.frame()


g2<-ggplot() +
  geom_point(data=res.lat,aes(x=median.temp,y=median.lat)) +
  geom_line(data=res.lat,aes(x=median.temp,y=median.lat)) +
  geom_line(data=res.lat,aes(x=median.temp,y=q1),linetype=3) +
  geom_line(data=res.lat,aes(x=median.temp,y=q3),linetype=3) +
  geom_ribbon(data=res.lat,aes(x=median.temp,ymin=q1,ymax=q3),alpha=0.3) +
  theme_bw() +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank()) +
  ylab("Abs. Latitude")


median.temp.dist<-NULL
max.temp.dist<-NULL
for (i in unique(res.env$bin)){
  max.temp.dist<-c(max.temp.dist,max(c(dist(res.env$Temperature[res.env$bin==i]))))
  median.temp.dist<-c(median.temp.dist,median(c(dist(res.env$Temperature[res.env$bin==i]))))
}
median(max.temp.dist)
median(median.temp.dist)

plot(res.env$median.temp[match(unique(res.env$bin),res.env$bin)],median.temp.dist,pch=2,ylim=c(0,15))
points(res.env$median.temp[match(unique(res.env$bin),res.env$bin)],max.temp.dist)

tp.df<-data.frame(temp=res.env$median.temp[match(unique(res.env$bin),res.env$bin)],median.t=median.temp.dist,max.t=max.temp.dist)
gextra<-ggplot() +
  geom_line(data=tp.df,aes(x=temp,y=median.t)) +
  geom_point(data=tp.df,aes(x=temp,y=median.t),shape=17) +
  geom_line(data=tp.df,aes(x=temp,y=max.t)) +
  geom_point(data=tp.df,aes(x=temp,y=max.t),shape=15) +
  theme_bw() +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank()) +
  ylab("Temperature difference")
gextra

res.ocean<-as.data.frame(table(res.env$bin,res.env$Ocean.region))
res.ocean$Var1<-as.integer(as.character(res.ocean$Var1))
res.ocean<-res.ocean %>%
  dplyr::rename(bin=Var1,Ocean=Var2,prop=Freq) %>%
  left_join(res.lat,by="bin")
res.ocean$Ocean<-gsub("\\[","",gsub("\\]","",sapply(strsplit(as.character(res.ocean$Ocean)," "),"[[",1)))
g3<-ggplot(data=res.ocean,aes(x=median.temp,y=prop,fill=Ocean)) +
  geom_area() +
  theme_bw() +
  scale_fill_manual(values=sushi.palette()) +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank()) +
  ylab("Number of samples")

# Complete plot --------------------------------------------------------------------------

g<-g3 / g2 / g1 + plot_layout(heights = c(2,2,6))
gB<-g3 / gextra / g1 + plot_layout(heights = c(2,2,6))
ggsave("../results/Fig4B_Varpart_EPI.pdf",g,height=8,width=5)
ggsave("../results/Fig4B_Varpart_EPI_B.pdf",gB,height=8,width=5)

# Plot with sampling hours on the panel --------------------------------------------------

res.hour<-res.env %>%
  group_by(bin,median.temp) %>%
  summarise(median.hour=median(hour),q1=quantile(hour,probs = 0.25),q3=quantile(hour,probs = 0.75)) %>%
  as.data.frame()

g2b<-ggplot() +
  geom_point(data=res.hour,aes(x=median.temp,y=median.hour)) +
  geom_line(data=res.hour,aes(x=median.temp,y=median.hour)) +
  geom_line(data=res.hour,aes(x=median.temp,y=q1),linetype=3) +
  geom_line(data=res.hour,aes(x=median.temp,y=q3),linetype=3) +
  geom_ribbon(data=res.hour,aes(x=median.temp,ymin=q1,ymax=q3),alpha=0.3) +
  theme_bw() +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank()) +
  ylab("Hour")
g<-g3 / g2b / g1 + plot_layout(heights = c(2,2,6))
ggsave("../results/Fig4B_Varpart_EPI_hour_in_panelB.pdf",g,height=8,width=5)

# Plot without interactions --------------------------------------------------------------

# Without interaction
g1<-ggplot() +
  geom_point(data=res.med %>% filter(Component.x!="interaction"),aes(x=median.temp,y=median,col=Component.x,shape=sign)) +
  #geom_line(data=res.med %>% filter(Component.x!="interaction"),aes(x=median.temp,y=median,col=Component.x)) +
  geom_smooth(data=res.med %>% filter(Component.x!="interaction"),aes(x=median.temp,y=median,col=Component.x),method="loess",span=0.1,se=F) +
  geom_ribbon(data=res.med %>% filter(Component.x!="interaction"),aes(x=median.temp,ymin=q1,ymax=q3,fill=Component.x),alpha=0.3) +
  #geom_text(data=res.med %>% filter(Component.x!="interaction"),aes(x=median.temp,y=3.2,label=sign),size=4) +
  theme_bw() +
  ylab("Squared euclidean distance") +
  xlab("Bins' median temperature") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c("#eaad53","#90e25a","gray")) +
  scale_fill_manual(values = c("#eaad53","#90e25a","gray")) +
  scale_shape_manual(values = c(3,19))

g<-g3 / g2 / g1 + plot_layout(heights = c(2,2,6))
ggsave("../results/Fig4B_Varpart_EPI_nointer.pdf",g,height=8,width=5)

# MES - Variance partitioning: Sample bins through temperature (each layers) -------------

# MES
res.mes<-varpart.sqr.euc.all(metaT.norm.match.log2[env.mat.match$Layer %in% c("MES"),],metaG.norm.match.log2[env.mat.match$Layer %in% c("MES"),],ratio.mat[env.mat.match$Layer %in% c("MES"),])
res.mes.comp.norm<-res.mes$components %>% gather("Component","value",-sample1,-sample2)
res.mes.comp.norm$Layer<-"MES"

mat.env<-env.mat.match[env.mat.match$Layer %in% c("MES"),]
n<-15
ordre<-order(mat.env$Temperature)
res<-NULL
for (i in 1:(nrow(mat.env)-n)){
  smpls<-mat.env$Barcode[ordre[i:(i+n-1)]]
  pos<-which(res.mes.comp.norm$sample1 %in% smpls & res.mes.comp.norm$sample2 %in% smpls)
  tmp<-res.mes.comp.norm[pos,]
  tmp$bin<-i
  tmp$median.temp<-median(mat.env$Temperature[ordre[i:(i+n-1)]])
  res<-rbind(res,tmp)
}


res.med<-res %>%
  group_by(bin,Component) %>%
  summarise(mean=mean(value),sd=sd(value),median=median(value),q1=quantile(value,probs = 0.25),q3=quantile(value,probs = 0.75)) %>%
  as.data.frame() %>%
  left_join(res,by="bin") %>%
  dplyr::select(-Component.y)

fun<-function(x){
  test<-wilcox.test(res$value[res$bin==x & res$Component=="Abundance"],res$value[res$bin==x & res$Component=="Expression"])
  test$p.value
}
pval.df<-data.frame(bin=unique(res.med$bin),p.val=sapply(unique(res.med$bin),fun))
pval.df$p.val.fdr<-p.adjust(pval.df$p.val,method = "bonferroni")
pval.df$sign<-sapply(pval.df$p.val.fdr,function(x){if (x<0.01) "P<0.05" else "NS"})

res.med<-left_join(res.med,pval.df,by="bin") %>%
  dplyr::select(-sample1,-sample2,-value) %>%
  unique()
res.med$sign[res.med$Component.x=="interaction"]<-"P<0.05"

# Plot distance partition with interactions ----------------------------------------------

g1<-ggplot() +
  geom_point(data=res.med,aes(x=median.temp,y=median,col=Component.x,shape=sign)) +
  geom_line(data=res.med,aes(x=median.temp,y=median,col=Component.x)) +
  #geom_smooth(data=res.med,aes(x=median.temp,y=median,col=Component.x),method="loess",span=0.1,se=F) +
  geom_ribbon(data=res.med,aes(x=median.temp,ymin=q1,ymax=q3,fill=Component.x),alpha=0.3) +
  #geom_text(data=res.med,aes(x=median.temp,y=3.2,label=sign),size=4) +
  theme_bw() +
  ylab("Component value") +
  xlab("Bins' median temperature") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c("#eaad53","#90e25a","gray")) +
  scale_fill_manual(values = c("#eaad53","#90e25a","gray")) +
  scale_shape_manual(values = c(3,19))

g1.later2<-g1 +
  scale_y_continuous(limits = c(-4,4)) +
  scale_shape_manual(values = c(19)) +
  theme(legend.position = "none") +
  labs(title="Mesopelagic")

# Plot distance partition with all panels ------------------------------------------------

res.env<-res %>%
  dplyr::select(sample1,sample2,bin,median.temp) %>%
  gather(key = "a",value="Barcode",-bin,-median.temp) %>%
  dplyr::select(-a) %>%
  left_join(env.mat.match,by="Barcode") %>%
  unique()

res.hour<-res.env %>%
  group_by(bin,median.temp) %>%
  summarise(median.lat=median(abs(Latitude)),q1=quantile(abs(Latitude),probs = 0.25),q3=quantile(abs(Latitude),probs = 0.75)) %>%
  as.data.frame()
g2<-ggplot() +
  geom_point(data=res.lat,aes(x=median.temp,y=median.lat)) +
  geom_line(data=res.lat,aes(x=median.temp,y=median.lat)) +
  geom_line(data=res.lat,aes(x=median.temp,y=q1),linetype=3) +
  geom_line(data=res.lat,aes(x=median.temp,y=q3),linetype=3) +
  geom_ribbon(data=res.lat,aes(x=median.temp,ymin=q1,ymax=q3),alpha=0.3) +
  theme_bw() +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank()) +
  ylab("Abs. Latitude")

res.ocean<-as.data.frame(table(res.env$bin,res.env$Ocean.region))
res.ocean$Var1<-as.integer(as.character(res.ocean$Var1))
res.ocean<-res.ocean %>%
  dplyr::rename(bin=Var1,Ocean=Var2,prop=Freq) %>%
  left_join(res.lat,by="bin")
res.ocean$Ocean<-gsub("\\[","",gsub("\\]","",sapply(strsplit(as.character(res.ocean$Ocean)," "),"[[",1)))
g3<-ggplot(data=res.ocean,aes(x=median.temp,y=prop,fill=Ocean)) +
  geom_area() +
  theme_bw() +
  scale_fill_manual(values=sushi.palette()) +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank()) +
  ylab("Number of samples")


g<-g3 / g2 / g1 + plot_layout(heights = c(2,2,6))
ggsave("../results/Fig4B_Varpart_MES.pdf",g,height=8,width=5)

# Plot distance partition without interactions -------------------------------------------

# Without interaction
g1<-ggplot() +
  geom_point(data=res.med %>% filter(Component.x!="interaction"),aes(x=median.temp,y=median,col=Component.x,shape=sign)) +
  geom_line(data=res.med %>% filter(Component.x!="interaction"),aes(x=median.temp,y=median,col=Component.x)) +
  #geom_smooth(data=res.med %>% filter(Component.x!="interaction"),aes(x=median.temp,y=median,col=Component.x),method="loess",span=0.1,se=F) +
  geom_ribbon(data=res.med %>% filter(Component.x!="interaction"),aes(x=median.temp,ymin=q1,ymax=q3,fill=Component.x),alpha=0.3) +
  #geom_text(data=res.med %>% filter(Component.x!="interaction"),aes(x=median.temp,y=3.2,label=sign),size=4) +
  theme_bw() +
  ylab("Squared euclidean distance") +
  xlab("Bins' median temperature") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c("#eaad53","#90e25a","gray")) +
  scale_fill_manual(values = c("#eaad53","#90e25a","gray"))

g<-g3 / g2 / g1 + plot_layout(heights = c(2,2,6))
ggsave("../results/Fig4B_Varpart_MES_nointer.pdf",g,height=8,width=5)

# EPI and MES - Plot supplementary figure with interactions ------------------------------

# Merge epi and MES

(g1.later | g1.later2) + plot_layout(widths = c(3,1))
ggsave("../results/Fig4B_Varpart.pdf",height=8,width=9)
