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

#Test proportion of specific nodes
  #Make unitary BLs and ultrametric
Strees_clado <- lapply(Strees, function(TREE){TREE$edge.length <- rep(1,length(TREE$edge.length));return(TREE)}) %>% lapply(chronos)
lapply(Strees_clado, function(TRE){
  #cophenetic(TRE)["Brooksia rostrata","Ritteriella retracta"]>cophenetic(TRE)["Ritteriella retracta","Cyclosalpa polae"] & 
    cophenetic(TRE)["Cyclosalpa polae","Cyclosalpa sewelli"]<cophenetic(TRE)["Cyclosalpa polae","Cyclosalpa quadriluminis"]
    #cophenetic(TRE)["Salpa aspera","Cyclosalpa polae"]<cophenetic(TRE)["Ritteriella retracta","Cyclosalpa polae"] &
    }) %>% unlist() %>% table() -> cr; print(cr/30.01)

par(mfrow=c(1,1),mar=c(0,0,0,0), oma=c(0,0,0,0))
#Append species to phylogeny based on taxonomy and morphological similarity -- ARBITRARY EDGE LENGTHS
ori_phylo <- chronos(tree_salp); ori_phylo$tip.label <- (str_replace_all(tree_salp$tip.label, " ", "_"))
extended_phylo <- add.species.to.genus(ori_phylo, "Pegea sp.", "Pegea")
extended_phylo <-  bind.tip(extended_phylo, "Helicosalpa younti", position=0.15, where = 7, edge.length = 0.15+0.099358)
extended_phylo <-  bind.tip(extended_phylo, "Cyclosalpa pinnata", position=0.02, where = 36, edge.length = 0.02712)
extended_phylo <-  bind.tip(extended_phylo, "Pegea socia", position=0.02, where = 25, edge.length = 0.02)

plot(extended_phylo)
edgelabels(extended_phylo$edge.length)
plot(extended_phylo); nodelabels(); tiplabels()

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

dvsz_fiji <- read.csv("DV_Zooid.stolon.angle.tsv", sep="\t", stringsAsFactors = F)[,c(5,8,12)]
names(dvsz_fiji)[3] <- "DVZS_Angle"

PEGEA_tree = extended_phylo
PEGEA_tree$tip.label[which(PEGEA_tree$tip.label == "Pegea confoederata")] <- "Pegea socia"
pruned_dvsz <- dvsz_fiji[which(dvsz_fiji$Species %in% PEGEA_tree$tip.label),]
dvsz_tree <- drop.tip(PEGEA_tree, which(!(PEGEA_tree$tip.label %in% unique(dvsz_fiji$Species))))

#### [0] #### HYDRODYNAMICS

#ONLY Taxa with Speed

kineTraits <- unique_traits[!is.na(unique_traits$Speed..cm.s.),]
rownames(kineTraits) <- kineTraits$Species

kineTree <- drop.tip(tree_salp_pruned, 
         which(tree_salp_pruned$tip.label %in% 
                 traits$Species[is.na(traits$Speed..cm.s.)]))
contColors <- function(trait, color1, color2){
  C <- kineTraits[,trait]
  cM <- contMap(kineTree, setNames(C, kineTraits$Species))
  cM$cols[1:length(cM$cols)]<-colorRampPalette(c(color1,color2), space="Lab")(length(cM$cols))
  plot(cM, legend=FALSE)
  add.color.bar(0.5, cM$cols,title=trait,prompt=FALSE,x=0, y=2-0.08*(Ntip(cM$tree)-1),lwd=4,fsize=0.8)
}

contColors(names(kineTraits)[15], "white", "black") #speed cm
contColors(names(kineTraits)[18], "white", "black") #speed z/s
contColors(names(kineTraits)[19], "white", "black") #speed colony

phenogram(kineTree, setNames(kineTraits$Speed..colonysize.min.,kineTraits$Species))
phenogram(kineTree, setNames(kineTraits$Speed..zooid.s,kineTraits$Species))
phenogram(tree_salp_pruned, setNames(unique_traits$DV.Zooid..Stolon.angle,unique_traits$Species))

#Front drag
kineTraits[,c(2,8,13)] %>% mutate(baseCSA = Zooid.area.cross.section.packing/Zooid.number)->CSA
CSA = mutate(CSA, Species = rownames(CSA))
CSAcum=CSA
for(r in 1:nrow(CSA)){
  species <- CSA[r,]
  for(z in 3:60){
    ZN <- species
    ZN$Zooid.number <- z
    if(ZN)
    ZN$Zooid.area.cross.section.packing <- z*ZN$baseCSA
    CSAcum <- rbind(CSAcum,ZN)
  }
}

CSAspeed <- left_join(CSA, kineTraits[,c(1,15,18,19)])

#ggplot(CSAspeed %>% group_by(Architecture), aes(Zooid.number, Zooid.area.cross.section.packing*baseCSA))+geom_point(aes(color=Architecture), alpha=0.5)+theme_bw()
ggplot(CSAspeed, aes(Architecture, Zooid.area.cross.section.packing)) +
  geom_violin(aes(color=Architecture,fill=Architecture)) +
  geom_text(label=CSAspeed$Species, vjust = 1)+ 
  scale_color_manual(values=setNames(c("orange","chartreuse","green", "cornflowerblue", "turquoise", "turquoise", "cyan"), unique(CSAspeed$Architecture)))+
  scale_fill_manual(values=setNames(c("orange","chartreuse","green", "cornflowerblue", "turquoise", "turquoise", "cyan"), unique(CSAspeed$Architecture)))+
  theme_bw()

ggplot(CSAspeed, aes(Zooid.area.cross.section.packing, Speed..colonysize.min.)) +
  geom_point(aes(color=Architecture)) +
  geom_text(label=CSAspeed$Species, vjust = 0.9, hjust=-0.1)+ 
  scale_color_manual(values=setNames(c("orange","chartreuse","green", "cornflowerblue", "turquoise", "turquoise", "cyan"), unique(CSAspeed$Architecture)))+
  theme_bw()

ggplot(kineTraits, aes(DV.Zooid..Stolon.angle, SN.Jet.Motion.angle)) +
  geom_point(aes(color=Architecture)) +
  geom_text(label=CSAspeed$Species, vjust = 0.9, hjust=-0.1)+ 
  scale_color_manual(values=setNames(c("orange","chartreuse","green", "cornflowerblue", "turquoise", "turquoise", "cyan"), unique(CSAspeed$Architecture)))+
  theme_bw()

#Power scaling
power <- full_join(CSA, kineTraits)
power %>% mutate(JetPower = Zooid.number.kine*Zooid.length.kine..cm.) %>% 
ggplot(aes(JetPower, Speed..cm.s.)) +
  geom_point(aes(color=Architecture)) +
  geom_text(label=power$Species, vjust = 0.9, hjust=-0.1)+ 
  theme_bw()

glm(Speed..colonysize.min. ~ JetPower, data=power %>% mutate(JetPower = Zooid.number.kine*Zooid.length.kine..cm.)) %>% 
  summary()

cbind(power, glm(Speed..colonysize.min. ~ JetPower, data=power %>% mutate(JetPower = Zooid.number.kine*Zooid.length.kine..cm.)) %>% 
        .$residuals)->power
names(power)[31] <- "jet_scaling"

ggplot(power, aes(Architecture, jet_scaling)) +
  geom_violin(aes(color=Architecture, fill=Architecture)) +
  geom_text(label=power$Species)+ 
  scale_color_manual(values=setNames(c("orange","chartreuse","green", "cornflowerblue", "turquoise", "turquoise", "cyan"), unique(power$Architecture)))+
  theme_bw()+
  scale_fill_manual(values=setNames(c("orange","chartreuse","green", "cornflowerblue", "turquoise", "turquoise", "cyan"), unique(power$Architecture)))
  

ggplot(kineTraits[-15,], aes(x=Zooid.number*log((Zooid.length..mm.^3)), y = Speed..colonysize.min.)) +
  geom_line(aes(color=Architecture), cex=2) +
  theme_bw()+
  scale_color_manual(values=setNames(c("turquoise","purple","magenta", "yellow", "orange", "green", "red"), unique(unique_traits$Architecture)))+
  theme(axis.text.x = element_text(angle = 90))

#Jet angle, colony speed - color is architecture
ggplot(kineTraits, aes(SN.Jet.Motion.angle, Speed..colonysize.min.)) +
  geom_point(aes(color=Architecture)) +
  geom_text(label=kineTraits$Species, vjust = 0.9, hjust=-0.1)+ 
  theme_bw()


#Zooid angle, colony speed - color is architecture
ggplot(kineTraits, aes(zDV.Zooid..Stolon.angle, Speed..colonysize.min.)) +
  geom_point(aes(color=Architecture)) +
  geom_text(label=kineTraits$Species, vjust = 0.9, hjust=-0.1)+
  xlim(c(0, 110))+
  theme_bw()

#CSA, colony speed - color is architecture
ggplot(kineTraits, aes(Zooid.area.cross.section.packing, Speed..colonysize.min.)) +
  geom_point(aes(color=Architecture)) +
  geom_text(label=kineTraits$Species, vjust = 0.9, hjust=-0.1)+
  xlim(c(0, 40))+
  theme_bw()

#packing, colony speed - color is architecture
ggplot(kineTraits, aes(Zooid.length.packing, Speed..colonysize.min.)) +
  geom_point(aes(color=Architecture)) +
  geom_text(label=kineTraits$Species, vjust = 0.9, hjust=-0.1)+
  xlim(c(1, 10))+
  theme_bw()

#COT, species - color is architecture
ggplot(kineTraits %>% filter(!is.na(Mean.Net.Cost..mgO2.ml.m.)), aes(x=Speed..colonysize.min., y=Mean.Net.Cost..mgO2.ml.m.)) +
  geom_point(aes(color=Species, size=Speed..cm.s.)) +
  theme_bw()

#COT vs %cost, species - color is architecture
ggplot(kineTraits %>% filter(!is.na(Mean.Net.Cost..mgO2.ml.m.)), aes(x=X..Cost.from.Swimming, y=Mean.Net.Cost..mgO2.ml.m.)) +
  geom_point(aes(color=Species)) +
  theme_bw()

ggplot(kineTraits %>% filter(!is.na(Mean.Net.Cost..mgO2.ml.m.)), aes(x=Mean.Gross.COT..mgO2.ml.m., y=DV.Zooid..Stolon.angle)) +
  geom_point(aes(color=Species)) +
  theme_bw()

PMScolor <- function(data, tree, x, y, chromo){
  data <- data[match(tree$tip.label, data$Species),]
  tip.cols <- brewer.pal(8,"Set1")[data[,chromo] %>%  as.character() %>% as.factor() %>% as.numeric()]
  tip.cols <- c()
  names(tip.cols)<-tree$tip.label
  cols<-c(tip.cols[tree$tip.label],rep("black",tree$Nnode))
  cols<-c(tip.cols[tree$tip.label],rep("black",tree$Nnode))
  names(cols)<-1:(length(tree$tip)+tree$Nnode)
  phylomorphospace(tree,as.matrix(data[,c(x,y)]),control=list(col.node=cols),label="horizontal", xlab=colnames(data)[x], ylab=colnames(data)[y])
}
par(oma=rep(4,4))
names(kineTraits)
simZDSV<-make.simmap(kineTree,setNames(kineTraits$Architecture,kineTraits$Species),nsim=25,model="ER")

phylomorphospace(kineTraits, simZDSV[[9]], 29, 19, 2)

#All taxa

ggplot(unique_traits %>% filter(!is.na(Mean.Gross.COT..mgO2.ml.m.)), aes(x=Species, y = Mean.Gross.COT..mgO2.ml.m.)) +
  geom_point(aes(color=Architecture)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))

ggplot(unique_traits[-17,], aes(x=Zooid.number, y = (Colony.cross.sectional.area...mm2./(Zooid.width..mm.^2)))) +
  geom_line(aes(color=Architecture), cex=2) +
  theme_bw()+
  scale_color_manual(values=setNames(c("turquoise","purple","magenta", "yellow", "orange", "green", "red"), unique(unique_traits$Architecture)))+
  theme(axis.text.x = element_text(angle = 90))



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
fitDiscrete(tree_salp_pruned, morph, model="meristic", symmetric=T)$opt$aicc %>% print()
fitDiscrete(tree_salp_morph, morph, model="meristic", symmetric=F)$opt$aicc %>% print()

  #SIMMAP plotting
# Morph <- c(setNames("Linear", "Ihlea punctata"),setNames("Oblique", "Thalia orientalis"),setNames("Linear","Salpa younti"),
#            setNames(unique_traits[,2],unique_traits$Species), setNames("Bipinnate","Ritteriella amboinensis"), 
#            setNames("Linear","Salpa thompsoni"), setNames("Linear","Metcalfina hexagona"), setNames("Whorl","Cyclosalpa bakeri"),
#            setNames("Helical","Helicosalpa virgula"))
# Mphylo <- drop.tip(extended_phylo, which(extended_phylo$tip.label=="Cyclosalpa floridiana" | extended_phylo$tip.label=="Brooksia lacromae"))
# Mphylo <- drop.tip(extended_phylo, which(!(extended_phylo$tip.label %in% names(Morph))))


simmorph<-make.simmap(tree_salp_morph,morph,nsim=100,model="ER")
par(ask=F)
obj_t<-summary(simmorph,plot=FALSE)
cols_t<-setNames(c("turquoise", "magenta","gold","orange","red","green","purple"),mapped.states(simmorph)[,1])
plot(obj_t,colors=cols_t,fsize=0.8,cex=c(0.9,0.5), ftype="i")
add.simmap.legend(colors=cols_t,x=0, y=8 ,prompt=FALSE,fsize=0.9,)

  #ACE plotting
fitMorph <- ace(morph,tree_salp_morph,"discrete")
plotTree(tree_salp_morph,fsize=0.8,ftype="i", offset=1)
cols=brewer.pal(8,"Set1")[1:7]
nodelabels(node=1:tree_salp_morph$Nnode+Ntip(tree_salp_morph),pie=fitMorph$lik.anc,piecol=cols,cex=0.7)
tiplabels(pie=to.matrix(morph,morph), piecol=cols,cex=0.5)

  #corHMM
#hidMorph <- corHMM(tree_salp_morph,cbind(names(morph),morph),rate.cat=1)
#hidMorph$solution[is.na(hidMorph$solution)] <- 0
#diag(hidMorph$solution) <- -rowSums(hidMorph$solution)

#MorphSimmap <- makeSimmap(tree_salp_morph,cbind(names(morph),morph), model=hidMorph$solution, rate.cat=1)
#plotSimmap(MorphSimmap[[1]], fsize = 0.5)

#### [2] #### Whole colony architecture EXCLUDING Helical for being autapomorphic #####
morph_nb <- morph[which(morph!="Helical")]
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
    fitDiscrete(tree_salp_morph, T_I, model=models[m])$opt$aicc %>% print()
    #fitDiscrete(tree_salp_morph, T_I, model=models[m]) %>% print()
    #simT<-make.simmap(tree_salp_morph,T_I,nsim=100,model=models[m],pi=table(T_I)/sum(table(T_I)))
    #par(ask=F)
    #obj_t<-summary(simT,plot=FALSE)
    #cols_t<-setNames(palette()[1:2],mapped.states(simT)[,1])
    #plot(obj_t,colors=cols_t,fsize=0.8,cex=c(0.9,0.5), ftype="i")
    #add.simmap.legend(colors=cols_t,x=0, y=4,prompt=FALSE,fsize=0.9)
    #ace(T_I, tree_salp_morph, type = "discrete", model=models[m]) %>% print()
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

#write.nexus.data(binTransits,"binDevelopmentalTransits.nex")
#write.table(binTransits,"binDevelopmentalTransits.tsv",sep="\t",row.names = F)
#write.tree(tree_salp_morph,"Chordata18Srev_morph.tree")

#MultiPhylo
for(ch in 2:ncol(binTransits)){
  print(names(binTransits)[ch])
  lapply(STree_morph, function(t){phylosig(t, setNames(binTransits[,ch], binTransits$Species))}) %>% unlist() %>% summary() %>% print()
  phylosig(tree_salp_morph, setNames(binTransits[,ch], binTransits$Species))%>% print()
  lapply(STree_morph, function(t){fastAnc(t, setNames(binTransits[,ch], binTransits$Species))["20"]}) %>% unlist() %>% summary() %>% print()
  fastAnc(tree_salp_morph, setNames(binTransits[,ch], binTransits$Species))["20"] %>% print()
}

lapply(STree_morph, function(t){phylosig(t, setNames(binTransits$WC, transits$Species))}) %>% unlist() %>% summary()
lapply(STree_morph, function(t){fastAnc(t, setNames(binTransits$OL, transits$Species))["20"]}) %>% unlist() %>% summary()

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
    #temperate
simLat<-make.simmap(tree_salp_lat,temperate,nsim=100)
par(ask=F)
obj<-summary(simLat,plot=FALSE)
cols<-setNames(palette()[1:3],mapped.states(simLat)[,1])
plot(obj,colors=cols,fsize=0.8,cex=c(0.9,0.5), ftype="i")
add.simmap.legend(colors=cols,x=0, y=4,prompt=FALSE,fsize=0.9,)

  #Tropical
simTro<-make.simmap(tree_salp_lat,tropical,nsim=100)
par(ask=F)
obj<-summary(simTro,plot=FALSE)
cols<-setNames(palette()[1:3],mapped.states(simLat)[,1])
plot(obj,colors=cols,fsize=0.8,cex=c(0.9,0.5), ftype="i")
add.simmap.legend(colors=cols,x=0, y=4,prompt=FALSE,fsize=0.9,)

  #Binary IS.TEMPERATE
temperate_bin <- setNames(as.numeric(as.factor(temperate)),names(temperate))-2
temperate_bin[temperate_bin==-1] <- 1
fastfitLat <- fastAnc(tree_salp_lat, temperate_bin)
  #logistic regression
temperate_binPrun = temperate_bin[which(names(temperate_bin)%in%names(dvsz))]
dvszPrun = dvsz[which(names(dvsz)%in%names(temperate_bin))]
treeIJ = drop.tip(tree_salp, which(!(tree_salp$tip.label %in% names(dvszPrun))))
datIJ = data.frame(temp = temperate_binPrun, angle = dvszPrun)
logit = glm(temp~angle,data=datIJ, family="binomial")
print(logit)
library(phylolm)
phylogit = phyloglm(temp~angle,phy=treeIJ,data=datIJ,boot=100, btol=30)
print(phylogit)

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
angleTree <- dvsz_tree
mean_dvzs = aggregate(DVZS_Angle ~ Species, data = pruned_dvsz, mean)
se_dvzs = aggregate(DVZS_Angle ~ Species, data = pruned_dvsz, function(x){sd(x)/sqrt(length(x))})

dvsz <- setNames(mean_dvzs$DVZS_Angle, mean_dvzs$Species)
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
  tree_i$tip.label[which(tree_i$tip.label == "Pegea confoederata")] <- "Pegea socia"
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

#%>% unlist() %>% summary()

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

  #OUWie Model comparison dvsz vs temperate
library(OUwie)
regimetree = treeIJ
tempAnc = ace(temperate_binPrun, regimetree, type = "discrete")$lik.anc
HypTempAnc = 1:nrow(tempAnc)
for(row in 1:nrow(tempAnc)){
  SCR = scale(tempAnc[row,])
  print(SCR)
  HypTempAnc[row] <- rownames(SCR)[which(SCR[,1] == max(SCR[,1]))] %>% print()
}
regimetree$node.label = as.factor(HypTempAnc)
#regimetree$node.label = as.factor(rep(1,length(HypTempAnc)))
plotTree(regimetree)
nodelabels(text=regimetree$node.label,frame="none",adj=c(1.6,-0.45), cex=0.6);tiplabels(text=temperate_binPrun, frame="none", cex=0.6, adj=c(0,2))
OUwie_matrix = as.data.frame(matrix(ncol=7, nrow=1))
names(OUwie_matrix) = c("Character", "Ntips", "dAICc_BM", "dAICc_BMS", "dAICc_OU1", "dAICc_OUm",  "dAICc_OUv")
library(magrittr)
Cdata = data.frame(cbind(rownames(datIJ), datIJ$temp, datIJ$angle))
names(Cdata) = c("species", "regime", "trait")
Cdata$trait = as.numeric(Cdata$trait)
Cdata$species = as.character(Cdata$species)
Cdata$regime = as.character(Cdata$regime)
BM <- OUwie(regimetree, Cdata, model="BM1")
BMs <- OUwie(regimetree, Cdata, model="BMS")
OU1 <- OUwie(regimetree, Cdata, model="OU1")
OUm <- OUwie(regimetree, Cdata, model="OUM")
OUmv <- OUwie(regimetree, Cdata, model="OUMV")
OUwie_matrix[1,1] = "dvsz"
OUwie_matrix[1,2] = nrow(Cdata)
OUwie_matrix[1,3] = BM$AICc - min(c(BM$AICc, BMs$AICc, OU1$AICc, OUm$AICc, OUmv$AICc))
OUwie_matrix[1,4] = BMs$AICc - min(c(BM$AICc, BMs$AICc, OU1$AICc, OUm$AICc, OUmv$AICc))
OUwie_matrix[1,5] = OU1$AICc - min(c(BM$AICc, BMs$AICc, OU1$AICc, OUm$AICc, OUmv$AICc))
OUwie_matrix[1,6] = OUm$AICc - min(c(BM$AICc, BMs$AICc, OU1$AICc, OUm$AICc, OUmv$AICc))
OUwie_matrix[1,7] = OUmv$AICc - min(c(BM$AICc, BMs$AICc, OU1$AICc, OUm$AICc, OUmv$AICc))
OUwie_matrix[,3:7] <- apply(OUwie_matrix[,3:7], 2, round, 3)

#SURFACE Convergence?
surfaceALL <- function(data,tree){
  Tree <- nameNodes(drop.tip(tree, which(!(tree$tip.label %in% rownames(data)))))
  olist <- convertTreeData(Tree, data)
  otree<-olist[[1]]
  odata<-olist[[2]]
  fwd<-surfaceForward(otree, odata, aic_threshold = 0, exclude = 0, verbose = FALSE, plotaic = FALSE)
  k<-length(fwd)
  fsum<-surfaceSummary(fwd)
  bwd<-surfaceBackward(otree, odata, starting_model = fwd[[k]], aic_threshold = 0, only_best = TRUE, verbose = FALSE, plotaic = FALSE)
  bsum<-surfaceSummary(bwd)
  kk<-length(bwd)
  print("N regimes BWD")
  bsum$n_regimes %>% return()
  par(mfrow=c(1,2))
  surfaceTreePlot(Tree, bwd[[kk]], labelshifts = T) %>% return()
  surfaceTraitPlot(data, bwd[[kk]], whattraits = c(1,2)) %>% return()
  newsim<-surfaceSimulate(Tree, type="hansen-fit", hansenfit=fwd[[k]]$fit, shifts=fwd[[k]]$savedshifts, sample_optima=TRUE)
  newout<-runSurface(Tree, newsim$dat, only_best = TRUE)
  newsum<-surfaceSummary(newout$bwd)
  newkk<-length(newout$bwd)
  print("N regimes SIM")
  newsum$n_regimes %>% return()
}

expected_angles[,1:6] <- sapply(expected_angles[,1:6],as.numeric)
surfaceALL(expected_angles[,5:6],angleTree)

#### [9] #### REVBAYES COLONY ARCHITECTURE EVOLUTION ################
RBPCMtree <- read.nexus("RB_PCM/morph_output/ase_morph.tree")
morphAncStates <- read.table("RB_PCM/morph_output/morph.states",sep='\t',header = T,stringsAsFactors = F)[,c(-1:-2)]
morphAncStates[morphAncStates==TRUE] <- "T"
apply(morphAncStates, 2, function(n){round(100*table(n)/sum(table(n)),2)->tabz;return(setNames(as.vector(tabz),names(tabz)))}) -> nodeStates
par(mfrow=c(4,5),mar=c(1,1,1,1))

lapply(nodeStates, function(v){
  if(names(v)==c("C","W")){colors=c("yellow","magenta")}
  if(names(v)==c("C", "O", "T", "W")){colors=c("yellow","cyan","purple", "magenta")}
  if(names(v)==c("L")){colors=c("green")}
  if(names(v)==c("B","L","O")){colors=c("red","green","cyan")}
  if(names(v)==c("B","L","O","T")){colors=c("red","green","cyan","purple")} 
  if(names(v)==c("B","L","O","T", "W")){colors=c("red","green","cyan","purple", "magenta")}
  if(names(v)==c("O","T", "W")){colors=c("cyan","purple", "magenta")}
  pie(v, col=colors, labels=" ")})
