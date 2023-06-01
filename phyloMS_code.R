library(tidyverse)
library(ape)
library(phytools)
library(reshape2)
library(data.table)
library(geiger)
library(bayou)
library(RColorBrewer)
library(phylolm)

setwd("~/salp_ecomorphology/")

#LOAD consensus tree
tree_salp <- read.nexus("phylogeny/RevBayes/TIMETREE_GUC-Mm+N+Sanger2_MUSCLE_output/TimeTree_GUC-Mm+N+Sanger2_MUSCLE_mcmc_MAP.tre")
  #Remove duplicate and ingroup undescribed salp species

    ### ATTENTION: MAKE SURE THESE MAXIMIZE BRANCH LENGTHS  ####
tree_salp <- drop.tip(tree_salp, c("HQ015387.1_Pegea_confoederata","HQ015397.1_Cyclosalpa_quadriluminis","HQ015391.1_Cyclosalpa_affinis",
                                   "HQ015396.1_Cyclosalpa_polae","HQ015398.1_Cyclosalpa_sewelli","FM244866.1_Iasis_cylindrica","HQ015399.1_Iasis_cylindrica",
                                   "HQ015402.1_Iasis_cylindrica","HQ015401.1_Iasis_cylindrica","HQ015413.1_Thalia_democratica","HQ015414.1_Thalia_democratica",
                                   "D14366.1_Thalia_democratica","HQ015410.1_Ritteriella_retracta","23_D37-Rsp-B-1_Ritteriella_sp", "24_D37-Rret-OS-1_Ritteriella_retracta",
                                   "HQ015404.1_Brooksia_rostrata","HQ015408.1_Salpa_maxima","HQ015406.1_Salpa_thompsoni", "HQ015377.1_Salpidae_gen._nov._sp._nov._A", 
                                   "FM244865.1_Ihlea_racovitzai", "KR057222.1_Brooksia_lacromae", "9_Helicosalpa_virgula_MarcHughes_specimen2",
                                   "26_Cyclosalpa_quadriluminis_D39CquaOS1", "27_Cyclosalpa_polae_D38CpolB1", "16_Ritteriella_amboinensis_D28RambOS1",
                                   "3_Ihlea_punctata_non-spotted_D24MhexB1", "20_Ihlea_punctata_D31IpunOS1-2", "6_Iasis_cf_cylindrica_yellow-tail_D22SyouB1"))
  #Remove accession codes and _
tree_salp$tip.label <- str_remove_all(tree_salp$tip.label, ".+\\..+?_")
tree_salp$tip.label <- str_replace_all(tree_salp$tip.label, "_", " ")
tree_salp$tip.label <- str_remove_all(tree_salp$tip.label, "^\\d+ ")
tree_salp$tip.label <- str_remove_all(tree_salp$tip.label, " D\\d+.+$")
tree_salp$tip.label <- str_remove_all(tree_salp$tip.label, "D\\d+.+? ")
  #Correct spellings
tree_salp$tip.label[which(tree_salp$tip.label == "Cyclosalpa floridana")] <- "Cyclosalpa floridiana"
  #Remove outgroups
tree_salp <- drop.tip(tree_salp, c("Pyrosomella verticillata", "Pyrosoma atlanticum", "Pyrosoma godeauxi","Pyrostremma spinosum", "Clavelina meridionalis", "Pycnoclavella aff. detorta", "Ascidia ceratodes", "Perophora sagamiensis","Megalodicopia hians", "Chelyosoma siboja", "Ciona intestinalis", "Molgula manhattensis", "Oikopleura dioica","Halocynthia igaboja", "Echinorhinus cookei", "Myxine glutinosa", "Branchiostoma floridae", "Doliolum denticulatum", "Doliolum nationalis")) 
#write.tree(tree_salp, "GUCMmNSanger_TimeTree_salp18Sphylo.tre")

#Load phylogenetic uncertainty tree set (3001 trees from RevBayes)
Strees <- read.tree("phylogeny/RevBayes/TOPOLOGY_GUC-Mm+N+Sanger2_MUSCLE_output/GUC-Mm+N+Sanger2_MUSCLE_18S.trees")
  #remove duplicate species and undescribed ingroup
Strees <- lapply(Strees, drop.tip, c("HQ015387.1_Pegea_confoederata","HQ015397.1_Cyclosalpa_quadriluminis","HQ015391.1_Cyclosalpa_affinis",
                                     "HQ015396.1_Cyclosalpa_polae","HQ015398.1_Cyclosalpa_sewelli","FM244866.1_Iasis_cylindrica","HQ015399.1_Iasis_cylindrica",
                                     "HQ015402.1_Iasis_cylindrica","HQ015401.1_Iasis_cylindrica","HQ015413.1_Thalia_democratica","HQ015414.1_Thalia_democratica",
                                     "D14366.1_Thalia_democratica","HQ015410.1_Ritteriella_retracta","23_D37-Rsp-B-1_Ritteriella_sp", "24_D37-Rret-OS-1_Ritteriella_retracta",
                                     "HQ015404.1_Brooksia_rostrata","HQ015408.1_Salpa_maxima","HQ015406.1_Salpa_thompsoni", "HQ015377.1_Salpidae_gen._nov._sp._nov._A", 
                                     "FM244865.1_Ihlea_racovitzai", "KR057222.1_Brooksia_lacromae", "9_Helicosalpa_virgula_MarcHughes_specimen2",
                                     "26_Cyclosalpa_quadriluminis_D39CquaOS1", "27_Cyclosalpa_polae_D38CpolB1", "16_Ritteriella_amboinensis_D28RambOS1",
                                     "3_Ihlea_punctata_non-spotted_D24MhexB1", "20_Ihlea_punctata_D31IpunOS1-2", "6_Iasis_cf_cylindrica_yellow-tail_D22SyouB1"))
  #Remove straneous characters and Accession codes from tip labels
Strees <- lapply(Strees, function(t){t$tip.label %>% 
    str_remove_all(".+\\..+?_") %>% 
    str_replace_all("_", " ") %>% 
    str_remove_all("^\\d+ ") %>% 
    str_remove_all(" D\\d+.+$") %>% 
    str_remove_all("D\\d+.+? ") -> t$tip.label; return(t)})
#reroot trees
Strees <- lapply(Strees, function(t){reroot(t, which(t$tip.label == "Branchiostoma floridae"))}) 
  #Remove outgroups
Strees <- lapply(Strees, drop.tip, c("Pyrosomella verticillata", "Pyrosoma atlanticum", "Pyrosoma godeauxi","Pyrostremma spinosum", "Clavelina meridionalis", "Pycnoclavella aff. detorta", "Ascidia ceratodes", "Perophora sagamiensis","Megalodicopia hians", "Chelyosoma siboja", "Ciona intestinalis", "Molgula manhattensis", "Oikopleura dioica","Halocynthia igaboja", "Echinorhinus cookei", "Myxine glutinosa", "Branchiostoma floridae", "Doliolum denticulatum", "Doliolum nationalis"))
  #correct spelling
Strees <- lapply(Strees, function(t){t$tip.label[which(t$tip.label == "Cyclosalpa floridana")] <- "Cyclosalpa floridiana"; return(t)})
#Ladderize for homogeneity in tip order
Strees <- lapply(Strees, ladderize)
  #Make ultrametric
Strees <- lapply(Strees, chronos)

#Quantify each unique tree topology
ape::unique.multiPhylo(Strees, use.tip.label = F)->Strees_Unique
lapply(Strees_Unique, function(t){length(lapply(Strees, function(Tr){all.equal.phylo(t,Tr,use.edge.length = F, use.tip.label = F)}) %>% unlist() %>% .[which(.==TRUE)])}) %>% unlist() %>% as.numeric() -> BSratios
BS_ratios <- (BSratios*100)/sum(BSratios)
round(BS_ratios,2) %>% sort() %>% .[length(BS_ratios):1]
Strees_Unique <- Strees_Unique[order(BS_ratios)[length(BS_ratios):1]]

#Plot each unique tree variant
par(mfrow=c(3,3),mar=c(0,0,0,0), oma=c(0,0,0,0))
lapply(Strees_Unique[1:9],function(t){plot.phylo(t, use.edge.length = F)})
#lapply(Strees_Unique,function(t){plot.phylo(t, use.edge.length = F, cex=0.2)})

par(mfrow=c(1,1),mar=c(0,0,0,0), oma=c(0,0,0,0))

#Load literature data
salplit <- read.csv("salplit.tsv", sep="\t", stringsAsFactors = F)
traits <- read.csv("SalpPreliminaryPass.tsv", sep="\t", stringsAsFactors = F)
  #correct spellings
salplit$Species[which(salplit$Species == "Iasis (Weelia) cylindrica")] <- "Iasis cylindrica"
salplit$Species[which(salplit$Species == "Soestia (Iasis) zonaria")] <- "Soestia zonaria"

  #prune tree
extended_phylo = tree_salp
#  write.tree(extended_phylo, "RB_PCM/GUCMmNSanger2_TimeTree_salp18Sphylo.tre")
#tree_salp_pruned <- drop.tip(tree_salp, which(!(tree_salp$tip.label %in% unique(traits$Species))))
tree_salp_pruned <- drop.tip(extended_phylo, which(!(extended_phylo$tip.label %in% unique(traits$Species))))
  #prune traits
pruned_traits <- traits[which(traits$Species %in% extended_phylo$tip.label),]
unique_traits <- unique(pruned_traits)

#### [1] #### Whole colony architecture #####
morph <- setNames(unique_traits$Architecture, unique_traits$Species)
morph <- c(morph,setNames("Linear", "Soestia zonaria"),
           setNames("Oblique", "Thalia orientalis"),
           setNames("Linear", "Salpa younti"),
           setNames("Linear", "Salpa thompsoni"),
           setNames("Bipinnate", "Ritteriella amboinensis"),
           setNames("Helical", "Helicosalpa younti"),
           setNames("Whorl","Cyclosalpa quadriluminis"))
tree_salp_morph <- drop.tip(extended_phylo, which(!(extended_phylo$tip.label %in% names(morph))))
STree_morph <- lapply(Strees, function(t){drop.tip(t,which(!(t$tip.label %in% names(morph))))})
morph <- morph[which(names(morph) %in% tree_salp_morph$tip.label)]
morph[match(tree_salp_morph$tip.label,names(morph))] -> morph

  #Model selection

models=c("ER","SYM","ARD")
for(i in 1:3){
  print(models[i])
  fitDiscrete(tree_salp_morph, morph, model=models[i])$opt$aicc %>% print() 
}

simmorph<-make.simmap(tree_salp_morph,morph,nsim=100,model="ER")
par(ask=F)
obj_t<-summary(simmorph,plot=FALSE)
cols_t<-setNames(c("turquoise", "magenta","gold","orange","red","green","purple"),mapped.states(simmorph)[,1])
plot(obj_t,colors=cols_t,fsize=0.8,cex=c(0.9,0.5), ftype="i")
add.simmap.legend(colors=cols_t,x=0, y=8 ,prompt=FALSE,fsize=0.9,)


#### [DVZS] ####  ANGLES ################
#DorsoVentral StolonZooid  angle

dvsz_fiji <- read.csv("DV_Zooid.stolon.angle.tsv", sep="\t", stringsAsFactors = F)[,c(5,8,12)]
names(dvsz_fiji)[3] <- "DVZS_Angle"

PEGEA_tree = extended_phylo
PEGEA_tree$tip.label[which(PEGEA_tree$tip.label == "Pegea confoederata")] <- "Pegea socia"
pruned_dvsz <- dvsz_fiji[which(dvsz_fiji$Species %in% PEGEA_tree$tip.label),]
dvsz_tree <- drop.tip(PEGEA_tree, which(!(PEGEA_tree$tip.label %in% unique(dvsz_fiji$Species))))
angleTree <- dvsz_tree
mean_dvzs = aggregate(DVZS_Angle ~ Species, data = pruned_dvsz, mean)
se_dvzs = aggregate(DVZS_Angle ~ Species, data = pruned_dvsz, function(x){sd(x)/sqrt(length(x))})
dvsz <- setNames(mean_dvzs$DVZS_Angle, mean_dvzs$Species) #named vector

# estimate ancestors
AncDVSZ<-fastAnc(angleTree,dvsz,CI=TRUE)
treePaint<-paintSubTree(angleTree,node=length(angleTree$tip)+1,"1")

# phenogram
trans<-as.character(floor(0:50/2))
trans[as.numeric(trans)<10]<- paste("0", trans[as.numeric(trans)<10],sep="")
for(i in 0:50){
  p<-i/length(trans)
  phenogram(angleTree,c(dvsz,(1-p)*AncDVSZ$CI95[,1]+p*AncDVSZ$ace), colors=setNames(paste("#0000ff",trans[i+1],sep=""),1), add=i>0)
  phenogram(angleTree,c(dvsz,(1-p)*AncDVSZ$CI95[,2]+p*AncDVSZ$ace), colors=setNames(paste("#0000ff",trans[i+1],sep=""),1), add=TRUE)
}
phenogram(angleTree,c(dvsz,AncDVSZ$ace),add=TRUE, colors=setNames("black",1))

phylosig(angleTree, dvsz, se=setNames(se_dvzs$DVZS_Angle, mean_dvzs$Species), test=TRUE)
fitContinuous(angleTree, dvsz, SE= setNames(se_dvzs$DVZS_Angle, mean_dvzs$Species))

#GEIGER Reversible Jump MCMC
#relaxed BM
rjmcmc.bm(angleTree, dvsz, ngen=10000, type="rbm", constrainSHIFT=1)
RBMdvsz <- load.rjmcmc("relaxedBM.result")
plot(x=RBMdvsz, par="shifts", burnin=0.25, edge.width=2)
RBMdvsz$log[which(RBMdvsz$log[,9]==max(RBMdvsz$log[,9])),8]
RBMdvsz$log[,9] %>% aicm()

rjmcmc.bm(angleTree, dvsz, ngen=10000, type="rbm", constrainSHIFT=2)
RBMdvsz <- load.rjmcmc("relaxedBM.result")
plot(x=RBMdvsz, par="shifts", burnin=0.25, edge.width=2)
RBMdvsz$log[which(RBMdvsz$log[,9]==max(RBMdvsz$log[,9])),8]
RBMdvsz$log[,9] %>% aicm()

rjmcmc.bm(angleTree, dvsz, ngen=10000, type="rbm", constrainSHIFT=3)
RBMdvsz <- load.rjmcmc("relaxedBM.result")
plot(x=RBMdvsz, par="shifts", burnin=0.25, edge.width=2)
RBMdvsz$log[which(RBMdvsz$log[,9]==max(RBMdvsz$log[,9])),8]
RBMdvsz$log[,9] %>% aicm()

rjmcmc.bm(angleTree, dvsz, ngen=10000, type="rbm", constrainSHIFT=4)
RBMdvsz <- load.rjmcmc("relaxedBM.result")
plot(x=RBMdvsz, par="shifts", burnin=0.25, edge.width=2)
RBMdvsz$log[which(RBMdvsz$log[,9]==max(RBMdvsz$log[,9])),8]
RBMdvsz$log[,9] %>% aicm() 

rjmcmc.bm(angleTree, dvsz, ngen=10000, type="rbm", constrainSHIFT=5)
RBMdvsz <- load.rjmcmc("relaxedBM.result")
plot(x=RBMdvsz, par="shifts", burnin=0.25, edge.width=2)
RBMdvsz$log[which(RBMdvsz$log[,9]==max(RBMdvsz$log[,9])),8]
RBMdvsz$log[,9] %>% aicm() #WInner

rjmcmc.bm(angleTree, dvsz, ngen=10000, type="rbm", constrainSHIFT=6)
RBMdvsz <- load.rjmcmc("relaxedBM.result")
plot(x=RBMdvsz, par="shifts", burnin=0.25, edge.width=2)
RBMdvsz$log[which(RBMdvsz$log[,9]==max(RBMdvsz$log[,9])),8]
RBMdvsz$log[,9] %>% aicm()

#BM + jumps
rjmcmc.bm(angleTree, dvsz, ngen=10000, type="jump-bm", constrainJUMP=4)
BMJdvsz <- load.rjmcmc("jump-BM.result")
plot(x=BMJdvsz, par="jumps", burnin=0.25, edge.width=2)
BMJdvsz$log[which(BMJdvsz$log[,9]==max(BMJdvsz$log[,9])),8]
BMJdvsz$log[,9] %>% aicm()

rjmcmc.bm(angleTree, dvsz, ngen=10000, type="jump-bm", constrainJUMP=3)
BMJdvsz <- load.rjmcmc("jump-BM.result")
plot(x=BMJdvsz, par="jumps", burnin=0.25, edge.width=2)
BMJdvsz$log[which(BMJdvsz$log[,9]==max(BMJdvsz$log[,9])),8]
BMJdvsz$log[,9] %>% aicm()

rjmcmc.bm(angleTree, dvsz, ngen=10000, type="jump-bm", constrainJUMP=2)
BMJdvsz <- load.rjmcmc("jump-BM.result")
plot(x=BMJdvsz, par="jumps", burnin=0.25, edge.width=2)
BMJdvsz$log[which(BMJdvsz$log[,9]==max(BMJdvsz$log[,9])),8]
BMJdvsz$log[,9] %>% aicm()  #WINNER!

rjmcmc.bm(angleTree, dvsz, ngen=10000, type="jump-bm", constrainJUMP=1)
BMJdvsz <- load.rjmcmc("jump-BM.result")
plot(x=BMJdvsz, par="jumps", burnin=0.25, edge.width=2)
BMJdvsz$log[which(BMJdvsz$log[,9]==max(BMJdvsz$log[,9])),8]
BMJdvsz$log[,9] %>% aicm()

pp.mcmc(phy=angleTree, d=dvsz)

#relaxed BM + jumps
rjmcmc.bm(angleTree, dvsz, ngen=10000, type="jump-rbm", constrainJUMP=2, constrainSHIFT=2)
RBMJdvsz <- load.rjmcmc("jump-relaxedBM.result")
plot(x=RBMJdvsz, par="jumps", burnin=0.25, edge.width=2)
dev.new()
plot(x=RBMJdvsz, par=c("shifts","jumps"), burnin=0.25, edge.width=2)
RBMJdvsz$log[which(RBMJdvsz$log[,9]==max(RBMJdvsz$log[,9])),8]
RBMJdvsz$log[,9] %>% aicm()

#MULTIPHYLO DV-St:Zoo
#Phylosig
Strees_PS = list()
for(i in 1:length(Strees_Unique)){
  tree_i = Strees_Unique[[i]]
  tree_i$tip.label[which(tree_i$tip.label == "Pegea confoederata")] <- "Pegea socia" #equivalent phylogenetic positions
  data_i <- dvsz_fiji[which(dvsz_fiji$Species %in% tree_i$tip.label),]
  tree_i <- drop.tip(tree_i, which(!(tree_i$tip.label %in% unique(data_i$Species))))
  mean_i = aggregate(DVZS_Angle ~ Species, data = data_i, mean)
  se_i = aggregate(DVZS_Angle ~ Species, data = data_i, function(x){sd(x)/sqrt(length(x))})
  mean_iv <- setNames(mean_i$DVZS_Angle, mean_i$Species)
  se_iv <- setNames(se_i$DVZS_Angle, se_i$Species)
  PS <- phylosig(tree_i, mean_iv, se=se_iv, test=TRUE)
  print(PS)
  Strees_PS[[i]] <- PS
}

###



