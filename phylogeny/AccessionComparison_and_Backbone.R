require(ape)


#Chordata Govindarajan Accession comparison

CH <- read.csv("phylogeny/chordata18S.fasta", header=F, stringsAsFactors = F)
CH[which(1:nrow(CH) %% 2 != 0),1] %>% str_remove_all("_.+") %>% sort() -> CHsort #52

GV <- read.csv("phylogeny/Govindarajan2011_18S.fasta", header=F, stringsAsFactors = F)
GV[which(1:nrow(GV) %% 2 != 0),1] %>% str_remove_all(" .+") %>% sort() -> GVsort #62

GVsort[which(GVsort %in% CHsort)] #48 of GV are in CH
CHsort[which(CHsort %in% GVsort)] #48 of CH are in GV
GVsort[which(!(GVsort %in% CHsort))] #14 of GV are NOT in CH
CHsort[which(!(CHsort %in% GVsort))] #4 of CH are NOT in GV

GVHtree <- read.newick("phylogeny/ML/Govindarajan2011_18S+Higaboja_aligned.fa.splits.nex")

#Backbone for Bayesian timetree
treeZ <- read.tree("phylogeny/ML/GovUNIONCho-AB921975Mman_aligned.fa.contree") #Load ML tree since RevBayes one has polytomies
treeX <- multi2di(treeZ) #remove a single polytomy somewhere - probably in the outgroups rooting
treeX$edge.length[1] <- 0.000001 #set arbitrary small branch length to artificial branch for polytomy
plot(treeX,cex=0.4);nodelabels(cex=0.4) #find rooting node (cephalochordata/vertebrata)
treeY <- reroot(treeX, 84) #set true root
treeU <- treeY %>% chronos() #make tree ultrametric
write.tree(treeU, "phylogeny/RevBayes/backbone_ML_Gov2011+Higa+Cint_MAFFT_rooted_multiD_ultra.tree")
