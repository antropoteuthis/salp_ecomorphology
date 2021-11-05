library(tidyverse)
library(ape)
library(phytools)
library(reshape2)

setwd("~/Documents/salp_ecomorphology/")
traits <- read.csv("salplit.tsv", sep="\t", stringsAsFactors = F)

ggplot(traits[which(traits$Variable=="Mean swimming speed cms"),], aes(x = Species, y = Value)) + geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(traits[which(traits$Variable=="Mean pulsation rate Hz"),], aes(x = Species, y = Value)) + geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(traits[which(traits$Variable=="Cost of locomotion J-kg-m"),], aes(x = Species, y = Value)) + geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#tree_salp <- read.newick("phylogeny/salps18Saligned.fa.treefile")
tree_salp <- read.newick("phylogeny/chordata18Saligned.fa.treefile")
tree_salp$tip.label <- str_remove_all(tree_salp$tip.label, ".+\\..+?_")
tree_salp$tip.label <- str_replace_all(tree_salp$tip.label, "_", " ")
tree_salp <- drop.tip(tree_salp, c(1:13, 45:52)) #drop outgroups
tree_salp <- drop.tip(tree_salp, c(29,28,24,23,21,20,19,16,12,13,9,7,4)) #drop duplicate species
tree_salp$tip.label[which(tree_salp$tip.label == "Iasis cylindrica")] <- "Iasis (Weelia) cylindrica"

#prune tree by data and make ultrametric
tree_salp_pruned <- drop.tip(tree_salp, which(!(tree_salp$tip.label %in% unique(traits$Species)))) %>% chronos()

pruned_traits <- traits[which(traits$Species %in% tree_salp_pruned$tip.label),]

cast_num <- dcast(pruned_traits[which(pruned_traits$Class=="number"),], Species~Variable, value.var="Value", fun.aggregate = function(x){mean(as.numeric(x), na.rm = T)})

#cast_cat <- dcast(traits[which(traits$Class=="category"),], Species~Variable, value.var="Value", fun.aggregate = unique)

morph <- setNames(pruned_traits[which(pruned_traits$Variable=="Chain architecture"),4], pruned_traits[which(pruned_traits$Variable=="Chain architecture"),1])
simmorph<-make.simmap(tree_salp_pruned,morph,model="ER",nsim=100)
par(ask=F)
obj<-summary(simmorph,plot=FALSE)
cols<-setNames(palette()[1:6],mapped.states(simmorph)[,1])
plot(obj,colors=cols,fsize=0.8,cex=c(0.9,0.5), ftype="i")
add.simmap.legend(colors=cols,x=0, y=4,prompt=FALSE,fsize=0.9,)

for(i in 2:ncol(cast_num)){
  CH_I=as.numeric(cast_num[,i])
  names(CH_I) = cast_num$Species
  CH_I= CH_I[!is.na(CH_I)]
  treeI = drop.tip(tree_salp_pruned,which(!(tree_salp_pruned$tip.label %in% names(CH_I))))
  class(treeI) = "phylo"
  PSIG <- phylosig(treeI, CH_I, test=T)
  print(names(cast_num)[i])
  PSIG$K %>% print()
  PSIG$P %>% print()
  length(CH_I) %>% print()
  #dotTree(treeI,CH_I)
  if(length(CH_I)>4){
    obj<-contMap(treeI,CH_I,plot=FALSE)
    grey<-setMap(obj,c("white","black"))
    par(mfrow=c(1,1))
    plot(grey,lwd=7,legend=2, leg.txt=names(cast_num)[i], fsize=c(1,0.5))
  }
}
