library(tidyverse)
library(ape)
library(phytools)
library(reshape2)
library(data.table)
library(geiger)
library(corHMM)
#library(diversitree)

setwd("~/Documents/salp_ecomorphology/")
traits <- read.csv("salplit.tsv", sep="\t", stringsAsFactors = F)
traits$Species[which(traits$Species == "Iasis (Weelia) cylindrica")] <- "Iasis cylindrica"
traits$Species[which(traits$Species == "Soestia (Iasis) zonaria")] <- "Soestia zonaria"

#tree_salp <- read.newick("phylogeny/salps18Saligned.fa.treefile")
tree_salp <- read.nexus("phylogeny/RevBayes/TIMETREE_Chordata_output/TimeTree_chordata_mcmc_MAP.tre")
tree_salp$tip.label <- str_remove_all(tree_salp$tip.label, ".+\\..+?_")
tree_salp$tip.label <- str_replace_all(tree_salp$tip.label, "_", " ")
tree_salp$tip.label[which(tree_salp$tip.label == "Cyclosalpa floridana")] <- "Cyclosalpa floridiana"
tree_salp <- drop.tip(tree_salp, c(32:52)) #drop outgroups
tree_salp <- drop.tip(tree_salp, c(1,2,6,9,14,16,18,19,20,22,23,27)) #drop duplicate species

#prune Literature data
tree_salp_pruned <- drop.tip(tree_salp, which(!(tree_salp$tip.label %in% unique(traits$Species))))
pruned_traits <- traits[which(traits$Species %in% tree_salp$tip.label),]
unique_traits <- unique(pruned_traits)

#### [1] #### Whole colony architecture #####
morph <- setNames(unique_traits[which(unique_traits$Variable=="Chain architecture"),4], unique_traits[which(unique_traits$Variable=="Chain architecture"),1])
morph <- c(morph,setNames("Linear", "Soestia zonaria"),setNames("Oblique", "Thalia orientalis"),setNames("Transversal", "Pegea confoederata"),setNames("Whorl","Cyclosalpa quadriluminis"))
tree_salp_morph <- drop.tip(tree_salp, which(!(tree_salp$tip.label %in% names(morph))))
morph <- morph[which(names(morph) %in% tree_salp_morph$tip.label)]
morph[match(tree_salp_morph$tip.label,names(morph))] -> morph
class(tree_salp_morph) <- "phylo"

  #Model selection

models=c("ER","SYM","ARD")
for(i in 1:3){
  print(models[i])
  fitDiscrete(tree_salp_morph, morph, model=models[i])$opt$aicc %>% print()
}
fitDiscrete(tree_salp_morph, morph, model="meristic", symmetric=T)$opt$aicc %>% print()
fitDiscrete(tree_salp_morph, morph, model="meristic", symmetric=F)$opt$aicc %>% print()

  #SIMMAP plotting
simmorph<-make.simmap(tree_salp_morph,morph,nsim=25,model="ARD",pi=table(morph)/sum(table(morph)))
par(ask=F)
obj_t<-summary(simmorph,plot=FALSE)
cols_t<-setNames(palette()[1:6],mapped.states(simmorph)[,1])
plot(obj_t,colors=cols_t,fsize=0.8,cex=c(0.9,0.5), ftype="i")
add.simmap.legend(colors=cols_t,x=0, y=4,prompt=FALSE,fsize=0.9,)

  #ACE plotting
fitMorph <- ace(morph,tree_salp_morph,"discrete")
plotTree(tree_salp_morph,fsize=0.8,ftype="i", offset=1)
cols=c("black","red","green","blue","cyan","pink")
nodelabels(node=1:tree_salp_morph$Nnode+Ntip(tree_salp_morph),pie=fitMorph$lik.anc,piecol=cols,cex=0.7)
tiplabels(pie=to.matrix(morph,sort(unique(morph))),piecol=cols,cex=0.5)
add.simmap.legend(colors=cols,prompt=FALSE,x=0,y=2,fsize=0.8)

  #corHMM
hidMorph <- corHMM(tree_salp_morph,cbind(names(morph),morph),rate.cat=1)
hidMorph$solution[is.na(hidMorph$solution)] <- 0
diag(hidMorph$solution) <- -rowSums(hidMorph$solution)

MorphSimmap <- makeSimmap(tree_salp_morph,cbind(names(morph),morph), model=hidMorph$solution, rate.cat=1)
plotSimmap(MorphSimmap[[1]], fsize = 0.5)

#### [2] #### Whole colony architecture EXCLUDING Bipinnate for being autapomorphic #####
morph_nb <- morph[which(morph!="Bipinnate")]
nbTree <- drop.tip(tree_salp, which(!(tree_salp$tip.label %in% names(morph_nb))))

  #Model selection
for(i in 1:3){
  print(models[i])
  fitDiscrete(nbTree, morph_nb, model=models[i])$opt$aicc %>% print()
}

fitDiscrete(nbTree, morph_nb, model="meristic", symmetric=T)$opt$aicc %>% print()
fitDiscrete(nbTree, morph_nb, model="meristic", symmetric=F)$opt$aicc %>% print()

nb_ordered <- fitMk(nbTree, x=morph_nb, model="ARD", ordered=T, pi=table(morph_nb)/sum(table(morph_nb)))
plot(nb_ordered,show.zeros=F)

  #SIMMAP plotting
simmorph_nb<-make.simmap(nbTree,morph_nb,nsim=25,model="ARD",pi=table(morph_nb)/sum(table(morph_nb)))
par(ask=F)
obj_nb<-summary(simmorph_nb,plot=FALSE)
cols_nb<-setNames(palette()[2:6],mapped.states(simmorph_nb)[,1])
plot(obj_nb,colors=cols_nb,fsize=0.8,cex=c(0.9,0.5), ftype="i")
add.simmap.legend(colors=cols_nb,x=0, y=4,prompt=FALSE,fsize=0.9,)

#### [3] #### BINARY Developmental transitions from Transversal budding ####
transits <- read.csv("Transitions_Salps.tsv", sep='\t', header = T, stringsAsFactors = F)[-4]
models=c("ER","ARD")
for(t in 2:5){
  T_I=transits[,t]
  names(T_I) = transits$Species
  for(m in 1:2){
    print(names(transits)[t])
    print(models[m])
    fitDiscrete(tree_salp_morph, T_I, model=models[m]) %>% print()
    simT<-make.simmap(tree_salp_morph,T_I,nsim=100,model=models[m],pi=table(T_I)/sum(table(T_I)))
    par(ask=F)
    obj_t<-summary(simT,plot=FALSE)
    cols_t<-setNames(palette()[1:2],mapped.states(simT)[,1])
    plot(obj_t,colors=cols_t,fsize=0.8,cex=c(0.9,0.5), ftype="i")
    add.simmap.legend(colors=cols_t,x=0, y=4,prompt=FALSE,fsize=0.9)
    ace(T_I, tree_salp_morph, type = "discrete", model=models[m]) %>% print()
  }
}

  #Discrete correlations between transitions - Bi-variate Mk model
  #don't work bc inapplicable combinations are allowed

par(mfrow=c(4,4))
for(t in 2:5){
  T_I=transits[,t]
  names(T_I) = transits$Species
  for(j in 2:5){
    if(t != j){
    T_J=transits[,j]
    names(T_J) = transits$Species
    fitPagel(tree_salp_morph,T_I,T_J) %>% plot(lwd.by.rate=TRUE)
    }
  }
}

  #Binary transitions as 1/0
binTransits <- transits
binTransits$TO <- as.numeric(as.factor(binTransits$TO))-1
binTransits$OL <- as.numeric(as.factor(binTransits$OL))
binTransits$OL[which(binTransits$OL == 2)] <- 0
binTransits$TW <- as.numeric(as.factor(binTransits$TW))-1
binTransits$WC <- as.numeric(as.factor(binTransits$WC))
binTransits$WC[which(binTransits$WC == 2)] <- 0

  #corHMM
hidTrans <- corHMM(tree_salp_morph,binTransits[,1:2],rate.cat=1)
hidTrans$solution[is.na(hidTrans$solution)] <- 0
diag(hidTrans$solution) <- -rowSums(hidTrans$solution)

HidSimmap <- makeSimmap(tree_salp_morph,binTransits[,1:2], model=hidTrans$solution, rate.cat=1)
plotSimmap(HidSimmap[[10]], fsize = 0.5)

corPAINT(tree_salp_morph, binTransits[,c(1:3)])

RDtransits <- rayDISC(tree_salp_morph, binTransits[,c(1,3)])

plotRECON(tree_salp_morph, RDtransits$states, title="OL")

ancRECON(tree_salp_morph, binTransits[,c(1:3)], p=param, rate.cat=NULL, ntraits=2, model="ARD")


#### [4] #### Colony morphologies as transition characters with exclusionary states (1/0/NA) ####
exc_transits <- read.csv("Exclusive_states.tsv", sep='\t', header = T, stringsAsFactors = F)
exc_transits <- as.data.frame(exc_transits %>% sapply(as.character))

  #model selection
models=c("ER","SYM","ARD")
par(mfrow=c(3,1))
for(t in c(2,3,5,6)){
  T_I=exc_transits[,t]
  names(T_I) = exc_transits$Species
  for(m in 1:length(models)){
    print(names(exc_transits)[t])
    print(models[m])
    fitDiscrete(tree_salp_morph, T_I, model=models[m]) %>% print()
    simT<-make.simmap(tree_salp_morph,T_I,nsim=100,model=models[m],pi=table(T_I)/sum(table(T_I)))
    par(ask=F)
    obj_t<-summary(simT,plot=FALSE)
    cols_t<-setNames(palette()[1:length(unique(T_I))],mapped.states(simT)[,1])
    plot(obj_t,colors=cols_t,fsize=0.8,cex=c(0.9,0.5), ftype="i")
    add.simmap.legend(colors=cols_t,x=0, y=4,prompt=FALSE,fsize=0.9)
    ace(T_I, tree_salp_morph, type = "discrete", model=models[m]) %>% print()
  }
}

  #model selection including meristics
models=c("ER","SYM","ARD","meristic","meristic")
for(t in c(2:6)){
  T_I=exc_transits[,t]
  names(T_I) = exc_transits$Species
  for(m in 1:length(models)){
    print(names(exc_transits)[t])
    print(models[m])
    if(m>3){
      fitDiscrete(tree_salp_morph, T_I, model=models[m],symmetric=ifelse(m==4,T,F))->fiti
    }
    else fitDiscrete(tree_salp_morph, T_I, model=models[m])->fiti
    fiti$opt$aicc %>% print()
    fiti$opt$k %>% print()
  }
}

  #SIMMAP plotting
par(mfrow=c(1,1))
for(t in c(2:6)){
  T_I=exc_transits[,t]
  names(T_I) = exc_transits$Species
    print(names(exc_transits)[t])
    simT<-make.simmap(tree_salp_morph,T_I,nsim=100,model="ER",pi=table(T_I)/sum(table(T_I)))
    par(ask=F)
    obj_t<-summary(simT,plot=FALSE)
    cols_t<-setNames(palette()[c(3,1,2)],mapped.states(simT)[,1])
    plot(obj_t,colors=cols_t,fsize=0.8,cex=c(0.9,0.5), ftype="i")
    add.simmap.legend(colors=cols_t,x=0, y=4,prompt=FALSE,fsize=0.9)
}


#### [5] #### #Developmental Transition pathways as multistate characters #####

transit_paths <- read.csv("Transitions_Paths.tsv", sep='\t', header = T, stringsAsFactors = F)

  #BIPPINATE PATH
bipQ <- rbind(c(0.1,0.1,0,0),c(0.1,0.1,0.1,0),c(0,0.1,0.1,0.1),c(0,0,0.1,0.1)) #Q matrix
rownames(bipQ)=colnames(bipQ)<-c("Bipinnate","Linear","Oblique","Transversal")
bipT <- transit_paths$Bipinnate_path
names(bipT) = transit_paths$Species

   #SIMMAP plotting
simBip<-make.simmap(tree_salp_morph,bipT,nsim=10,pi=table(bipT)/sum(table(bipT)),model="ARD")
par(ask=F)
obj_bip<-summary(simBip,plot=FALSE)
cols_bip<-setNames(palette()[1:length(unique(bipT))],mapped.states(simBip)[,1])
plot(obj_bip,colors=cols_bip,fsize=0.8,cex=c(0.9,0.5), ftype="i")
add.simmap.legend(colors=cols_bip,x=0, y=4,prompt=FALSE,fsize=0.9,)

   #Model selection
bip_meristic <- fitDiscrete(tree_salp_morph, bipT, model="meristic", symmetric=T)
models=c("ER","SYM","ARD")
for(i in 1:3){
  fitDiscrete(tree_salp_morph, bipT, model=models[i])->fiti
  fiti$opt$aicc %>% print()
  fiti$opt$k %>% print()
}
bip_ordered <- fitMk(tree_salp_morph, x=bipT,Q=bipQ, model="SYM", ordered=T, pi=table(bipT)/sum(table(bipT)))
plot(bip_ordered,show.zeros=F)
bip_ace <- ace(bipT, tree_salp_morph, type = "discrete", model="SYM")

  #Simplified Linear PATH (T-O-L)

linT <- transit_paths$Bipinnate_path
names(linT) = transit_paths$Species
linT[which(linT == "Linear")] <- "Oblique"

  #SIMMAP plotting
simLin<-make.simmap(tree_salp_morph,linT,nsim=25,pi=table(linT)/sum(table(linT)),model="ER")
par(ask=F)
obj_lin<-summary(simLin,plot=FALSE)
cols_lin<-setNames(palette()[1:length(unique(linT))],mapped.states(simLin)[,1])
plot(obj_lin,colors=cols_lin,fsize=0.8,cex=c(0.9,0.5), ftype="i")
add.simmap.legend(colors=cols_lin,x=0, y=4,prompt=FALSE,fsize=0.9,)

    #model selection
lin_meristic <- fitDiscrete(tree_salp_morph, linT, model="meristic", symmetric=T)
fitDiscrete(tree_salp_morph, linT, model="meristic", symmetric=F)
for(i in 1:3){
  fitDiscrete(tree_salp_morph, linT, model=models[i])->fiti
  fiti$opt$aicc %>% print()
  fiti$opt$k %>% print()
}

lin_ordered <- fitMk(tree_salp_morph, x=linT, model="SYM", ordered=T, pi=table(linT)/max(table(linT)))
plot(lin_ordered,show.zeros=F)

lin_ace <- ace(linT, tree_salp_morph, type = "discrete", model="SYM")

  #CLUSTER PATH
cluT <- transit_paths$Cluster_path
names(cluT) = transit_paths$Species
cluQ <- rbind(c(1,0,1),c(0,1,1),c(1,1,1))
rownames(cluQ)=colnames(cluQ)<-c("Cluster","Transversal","Whorl")

  #SIMMAP plotting

simClu<-make.simmap(tree_salp_morph,cluT,nsim=10,pi=table(cluT)/sum(table(cluT)),model="ER")
par(ask=F)
obj_clu<-summary(simClu,plot=FALSE)
cols_clu<-setNames(palette()[1:length(unique(cluT))],mapped.states(simClu)[,1])
plot(obj_clu,colors=cols_clu,fsize=0.8,cex=c(0.9,0.5), ftype="i")
add.simmap.legend(colors=cols_clu,x=0, y=4,prompt=FALSE,fsize=0.9,)

  #Model selection
fitDiscrete(tree_salp_morph, cluT, model="meristic", symmetric=T)
fitDiscrete(tree_salp_morph, cluT, model="meristic", symmetric=F)

models=c("ER","SYM","ARD")
for(i in 1:3){
  fitDiscrete(tree_salp_morph, cluT, model=models[i])->fiti
  fiti$opt$aicc %>% print()
  fiti$opt$k %>% print()
}

clu_ordered <- fitMk(tree_salp_morph, x=cluT, model="ER")
plot(clu_ordered,show.zeros=F)
clu_ace <- ace(cluT, tree_salp_morph, type = "discrete", model="SYM")

#### [6] #### LATITUDE RANGES Categorical #####

tropical <- setNames(unique_traits[which(unique_traits$Variable=="Tropical"),4], unique_traits[which(unique_traits$Variable=="Tropical"),1])
temperate <- setNames(unique_traits[which(unique_traits$Variable=="Temperate"),4], unique_traits[which(unique_traits$Variable=="Temperate"),1])
tree_salp_lat <- drop.tip(tree_salp, which(!(tree_salp$tip.label %in% names(temperate)))) %>% chronos()
temperate <- temperate[which(names(temperate) %in% tree_salp_lat$tip.label)]
temperate[match(names(temperate),tree_salp_lat$tip.label)] -> temperate

  #ACE plotting
fitLat <- ace(temperate,tree_salp_lat,"discrete")
plotTree(tree_salp_lat,fsize=0.8,ftype="i", offset=1)
cols=c("red","green","black")
nodelabels(node=1:tree_salp_lat$Nnode+Ntip(tree_salp_lat),pie=fitLat$lik.anc,piecol=cols,cex=0.7)
tiplabels(pie=to.matrix(tropical,sort(unique(tropical))),piecol=cols,cex=0.5)
add.simmap.legend(colors=cols,prompt=FALSE,x=0,y=2,fsize=0.8)

  #SIMMAP plotting
simLat<-make.simmap(tree_salp_lat,temperate,nsim=100)
par(ask=F)
obj<-summary(simLat,plot=FALSE)
cols<-setNames(palette()[1:3],mapped.states(simLat)[,1])
plot(obj,colors=cols,fsize=0.8,cex=c(0.9,0.5), ftype="i")
add.simmap.legend(colors=cols,x=0, y=4,prompt=FALSE,fsize=0.9,)

#### [7 #### Continuous characters from the literature #####
cast_num <- dcast(pruned_traits[which(pruned_traits$Class=="number"),], Species~Variable, value.var="Value", fun.aggregate = function(x){mean(as.numeric(x), na.rm = T)})

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
  }
}

#### [8] #### EXPECTED ANGLES ################

expected_angles <- read.csv("expected_angles.tsv",sep="\t",header = T, row.names = 1, stringsAsFactors = F)[,-c(7,8)] #remove variables with all zeroes

  #ContMaps phylogenetic signals
for(i in 2:ncol(expected_angles)){
  CH_I=as.numeric(expected_angles[,i])
  names(CH_I) = row.names(expected_angles)
  CH_I= CH_I[!is.na(CH_I)]
  treeI = drop.tip(tree_salp,which(!(tree_salp$tip.label %in% names(CH_I)))) %>% chronos()
  class(treeI) = "phylo"
  PSIG <- phylosig(treeI, CH_I, test=T)
  print(names(expected_angles)[i])
  PSIG$K %>% print()
  PSIG$P %>% print()
  length(CH_I) %>% print()
  #dotTree(treeI,CH_I)
  if(length(CH_I)>4){
    obj<-contMap(treeI,CH_I,plot=FALSE)
    grey<-setMap(obj,c("white","black"))
    par(mfrow=c(1,1))
    plot(grey,lwd=7,legend=2, leg.txt=names(expected_angles)[i], fsize=c(1,0.5))
  }
}

  #Model comparison on expected angles
AICdf = as.data.frame(matrix(ncol=6,nrow=2))
colnames(AICdf) = c("Variable", "white_noise", "starBM", "BM", "EB", "OU")
angle_tree <- chronos(tree_salp)
#angle_tree <- chronos(drop.tip(tree_salp, which(tree_salp$tip.label %in% c("Pegea confoederata","Pegea bicaudata","Cyclosalpa affinis","Cyclosalpa polae","Cyclosalpa sewelli","Cyclosalpa floridiana"))))
startree <- rescale(angle_tree, "lambda", 0)
for(c in 1:ncol(expected_angles)){
  C = expected_angles[which(rownames(expected_angles) %in% angle_tree$tip.label),c]
  names(C) = rownames(expected_angles)[which(rownames(expected_angles) %in% angle_tree$tip.label)]
  C = C[!is.na(C)]
  Ctree = drop.tip(angle_tree, which(!(angle_tree$tip.label %in% names(C))))
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
    model_matrix[j, "aicc"] <- temp_model$aicc
  }
  model_matrix[,"aicc_best"] <- min(model_matrix[,"aicc"])
  model_matrix[,"dAICc"] <- model_matrix[, "aicc"] - model_matrix[j, "aicc_best"]
  print(names(expected_angles)[c])
  string_c <- c(names(expected_angles)[c], model_matrix[,3])
  names(string_c) = colnames(AICdf)
  AICdf[c,] <- string_c
}
AICdf[,2:6] = apply(AICdf[,2:6], 2, as.numeric) %>% apply(2, function(x){round(x,3)})

#ace(C, Ctree, model="EB")
  #DorsoVentral StolonZooid Expected angle
DVSZ<-contMap(Ctree,C,plot=FALSE)
DVSZ<-setMap(DVSZ,c("white","black"))
par(mfrow=c(1,1))
plot(DVSZ,lwd=7,legend=2, leg.txt=names(expected_angles)[6], fsize=c(1,0.5))
phenogram(Ctree,C)

#SURFACE Convergence?
