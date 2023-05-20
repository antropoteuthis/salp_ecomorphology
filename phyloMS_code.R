library(tidyverse)
library(ape)
library(phytools)
library(reshape2)
library(data.table)
library(geiger)
#library(corHMM)
library(bayou)
library(surface)
library(RColorBrewer)
library(phylolm)
#library(diversitree)

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


#### [7] #### Continuous characters from the literature #####
#cast_num <- dcast(pruned_traits[which(pruned_traits$Class=="number"),], Species~Variable, value.var="Value", fun.aggregate = function(x){mean(as.numeric(x), na.rm = T)})
cast_num <- pruned_traits[,c(1,3:23,29)]
  #Phylogenetic signals and contMaps for each character
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
    phenogram(treeI, CH_I)
  }
}

#### [8] ####  ANGLES ################
#DorsoVentral StolonZooid Expected angle
angleTree <- drop.tip(tree_salp,which(!(tree_salp$tip.label %in% unique_traits$Species)))
dvsz <- setNames(unique_traits$DV.Zooid..Stolon.angle, unique_traits$Species)
#contMap
DVSZ<-contMap(angleTree, dvsz,plot=FALSE)
DVSZ<-setMap(DVSZ,c("white","orange","red","green"))
par(mfrow=c(1,1))
plot(DVSZ,lwd=7,legend=2, leg.txt="Dorsoventral Zooid-Stolon Angle", fsize=c(1,0.5))
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
RBMdvsz$log[,9] %>% aicm() #WInner

rjmcmc.bm(angleTree, dvsz, ngen=10000, type="rbm", constrainSHIFT=5)
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
BMJdvsz$log[,9] %>% aicm() #WINNER!

rjmcmc.bm(angleTree, dvsz, ngen=10000, type="jump-bm", constrainJUMP=2)
BMJdvsz <- load.rjmcmc("jump-BM.result")
plot(x=BMJdvsz, par="jumps", burnin=0.25, edge.width=2)
BMJdvsz$log[which(BMJdvsz$log[,9]==max(BMJdvsz$log[,9])),8]
BMJdvsz$log[,9] %>% aicm()

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
phylosig(angleTree, setNames(expected_angles$DVN_Stolon.Zooid, row.names(expected_angles)))
lapply(Strees, function(t){
  PS <- phylosig(t, setNames(expected_angles$DVN_Stolon.Zooid, row.names(expected_angles)))
  return(PS)
}) %>% unlist() %>% summary()

###

#Model comparison
AICdf = as.data.frame(matrix(ncol=6,nrow=2))
colnames(AICdf) = c("Variable", "white_noise", "starBM", "BM", "EB", "OU")
#angle_tree <- chronos(drop.tip(tree_salp, which(tree_salp$tip.label %in% c("Pegea confoederata","Pegea bicaudata","Cyclosalpa affinis","Cyclosalpa polae","Cyclosalpa sewelli","Cyclosalpa floridiana"))))
startree <- rescale(angleTree, "lambda", 0)
  C = dvsz
  Ctree = angleTree
  startree_C = drop.tip(startree, which(!(startree$tip.label %in% names(C))))
  model_matrix = matrix("NA", nrow = 5, ncol = 3)
  colnames(model_matrix) = c("aicc","aicc_best","dAICc")
  row.names(model_matrix) = c("white", "starBM", "BM", "EB", "OU")
  for(j in 1:dim(model_matrix)[1]){
    if(j==2){
      temp_model = fitContinuous(startree_C, C, model="BM")$opt
    }
    else{
      temp_model = fitContinuous(Ctree, C, model=row.names(model_matrix)[j])$opt
    }
    model_matrix = apply(model_matrix,2, as.numeric)
    row.names(model_matrix) = c("white", "starBM", "BM", "EB", "OU")
    model_matrix[j, "aicc"] <- temp_model$aic
  }
  model_matrix[,"aicc_best"] <- min(model_matrix[,"aicc"])
  model_matrix[,"dAICc"] <- model_matrix[, "aicc"] - model_matrix[j, "aicc_best"]
  print(names(expected_angles)[c])
  string_c <- c(names(expected_angles)[c], model_matrix[,3])
  names(string_c) = colnames(AICdf)
  AICdf[c,] <- string_c

AICdf[,2:6] = apply(AICdf[,2:6], 2, as.numeric) %>% apply(2, function(x){round(x,3)})

