library(data.table)
library(tidyverse)

# Load pairs
dat<-fread("../results/Coexpression_pairs_renamed.tsv",sep="\t",header=T,data.table = F,stringsAsFactors = F)

# Network
library(ggnetwork)
library(network)
dat<-dat[1:1000,]

nodes.df<-data.frame(lab=unique(c(dat$`Representative 1`,dat$`Representative 2`)))
nodes.df$node.id<-1:nrow(nodes.df)

edge.df<-data.frame(from=nodes.df$node.id[match(dat$`Representative 1`,nodes.df$lab)],to=nodes.df$node.id[match(dat$`Representative 2`,nodes.df$lab)],to_type=dat$`Representative 2 type`)
net<-network(edge.df,matrix.type = "edgelist", ignore.eval = FALSE)


plot.network(net,vertex.col=)

ggnet<-ggnetwork(net)
ggplot(ggnetwork(net), aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges() +
  geom_nodes(aes(color = to_type)) +
  theme_blank()
