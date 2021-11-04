library(tidyverse)
library(ape)
library(phytools)
library(reshape2)

setwd("~/Documents/salp_ecomorphology/")
traits <- read.csv("salplit.tsv", sep="\t", stringsAsFactors = F)

ggplot(traits[which(traits$Variable=="Mean swimming speed cms"),], aes(x = Species, y = Value)) + geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(traits[which(traits$Variable=="Pulsation rate Hz"),], aes(x = Species, y = Value)) + geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(traits[which(traits$Variable=="Cost of locomotion J-kg-m"),], aes(x = Species, y = Value)) + geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tree_salp <- read.newick("salpML.tre")
tree_salp <- rtree(length(unique(traits$Species))) ## BULLSHIT
tree_salp$tip.label <- unique(traits$Species) ## BULLSHIT

tree_salp_pruned <- drop.tip(tree_salp, which(!(tree_salp$tip.label %in% unique(traits$Species))))
pruned_traits <- traits[which(traits$Species %in% tree_salp_pruned$tip.label),]

cast_num <- dcast(traits[which(traits$Class=="number"),], Species~Variable, value.var="Value", fun.aggregate = function(x){mean(as.numeric(x), na.rm = T)})

cast_cat <- dcast(traits[which(traits$Class=="category"),], Species~Variable, value.var="Value", fun.aggregate = unique)
morph <- traits[which(traits$Variable=="Chain architecture"),c(1,4)]
names(morph)[2] <- "Chain architecture"
