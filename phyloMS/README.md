# 

 A new molecular phylogeny of salps (Tunicata: Thalicea: Salpida) and the evolutionary history of their colonial architecture
 DATA AND CODE FOR IOB MANUSCRIPT
---

Salps are marine pelagic urochordates with a complex life cycle including a solitary and colonial stage composed of asexually-budded individuals. These colonies develop into species-specific architectures with distinct zooid orientations, including transversal, oblique, linear, helical, and bipinnate chains; as well as whorls, and clusters. The evolutionary history of salp colony architecture has remained obscured due to the lack of a homology-based ontology to characterize architectures, as well as a lack of phylogenetic taxon sampling and resolution of critical nodes. We (1) collected and first-time sequenced eight species of salps, (2) inferred the phylogenetic relationships among salps, and (3) reconstructed the evolutionary history of salp colony architecture. We collected salp specimens via offshore SCUBA diving, dissected tissue samples, extracted their DNA, amplified their 18S gene, and sequenced them using Sanger technology. We inferred a new molecular phylogeny using both Maximum Likelihood and Bayesian approaches. Using this phylogeny, we reconstructed the ancestral states of colony architecture using a Bayesian ordered Markov model informed by the presence and absence of specific developmental mechanisms that lead to each architecture. We find that the ancestral salp architecture is either oblique or linear, with every other state being derived. Moreover, linear chains have evolved independently at least three times. While transversal chains are developmentally basal and hypothesized to be ancestral, our phylogenetic topology and reconstructions strongly indicate that they are evolutionarily derived through the loss of zooid torsion. These traits are likely critical to multijet locomotory performance and evolving under natural selection. Our work showcases the need to study the broader diversity of salp species in order to gain a comprehensive understanding of their organismal biology, evolutionary history, and ecological roles in pelagic ecosystems.

## Description of the data and file structure

### Tab-separated value spreadsheets (any spreadsheet software)

DV_Zooid.stolon.angle.tsv    Measurements of dorsoventral zooid-stolon angles including Filename, File.type, Source, Authorship, Species, Specimen_type, Specimen, Architecture, and Notes.
salplit.tsv    Data derived from the literature including colonial architecture as well as other ecological and physiological variables. 

### Data and logs from the RevBayes ancestral state estimation (text files)

morph.nex   Data input for ancestral state estimation derived from phyloMS_code.R pruning of salplit.tsv
GUCMmNSanger2_TimeTree_salp18Sphylo.tre   Ultrametric time tree input for ancestral state estimation derived from phyloMS_code.R pruning of TimeTree_GUC-Mm_N_Sanger2_MUSCLE_mcmc_MAP.tre
morph.states    Output from the analysis with multiple rates for each developmental transition using the mcmc_ase_freeK_RJ.Rev script
morph_single.states    Output from the analysis with a single rate parameter for all developmental transitions using the mcmc_ase_singleRate.Rev script
morph.log    Log from the analysis with multiple rates for each developmental transition using the mcmc_ase_freeK_RJ.Rev script
morph_single.log    Log from the analysis with a single rate parameter for all developmental transitions using the mcmc_ase_singleRate.Rev script
ase_morph.tree    Tree output from the analysis with multiple rates for each developmental transition using the mcmc_ase_freeK_RJ.Rev script
ase_morph_single.tree    Tree output from the analysis with a single rate parameter for all developmental transitions using the mcmc_ase_singleRate.Rev script

### Alignments (any alignment viewing software or text file software)

GUC-Mm_N_Sanger2_aligned.fa    MAFFT alignment output from: > mafft --auto --inputorder "Sanger/GUC-Mm+N+Sanger2.fa" > "GUC-Mm+N+Sanger2_aligned.fa"
GUC-Mm_N_Sanger2_aligned-Cbak.fa-gb   GBLOCKS postprocessing output of the MAFFT alignment from: > ./Gblocks  GUC-Mm+N+Sanger2_aligned.fa -t DNA -b5 h > GUC-Mm+N+Sanger2_aligned.fa-gb
GUC-Mm_N_Sanger2_aligned.fa-gb.uniqueseq.phy   GBLOCKS postprocessing Phylip output of the MAFFT alignment from: > ./Gblocks  GUC-Mm+N+Sanger2_aligned.fa -t DNA -b5 h > GUC-Mm+N+Sanger2_aligned.fa-gb

GUC-Mm_N_Sanger2_MUSCLE.afa    MUSCLE alignment output from: > ../muscle5.1.macos_arm64 -align Sanger/GUC-Mm+N+Sanger2.fa -output GUC-Mm+N+Sanger2_MUSCLE.afa
GUC-Mm_N_Sanger2-Cbak_MUSCLE.afa-gb    GBLOCKS postprocessing output of the MUSCLE alignment from: > ./Gblocks GUC-Mm+N+Sanger2_MUSCLE.afa -t DNA -b5 h > GUC-Mm+N+Sanger2_aMUSCLE.afa-gb
GUC-Mm_N_Sanger2_MUSCLE.afa-gb.uniqueseq.phy    GBLOCKS postprocessing Phylip output of the MUSCLE alignment from: > ./Gblocks GUC-Mm+N+Sanger2_MUSCLE.afa -t DNA -b5 h > GUC-Mm+N+Sanger2_aMUSCLE.afa-gb

### IQTree Maximum Likelihood outputs

Produced by: > iqtree -s ALIGNMENT -nt AUTO -bb 1000

Key to the IQTree outputs:

Files ended in .splits.nex are NEXUS representations of phylogenetic uncertainty
Files ended in .treefile contain the ML trees in NEWICK format (open in FigTree or iTOL) with bootstrap variants
Files ended in .contree contain the ML consensus tree in NEWICK format (open in FigTree or iTOL)
Files ended in .bionj contain the ML trees in BIONJ format
Files ended in .ckp.gz contain a checkpoint file in (gzip-compressed)
Files ended in .log contain the log file of the run
Files ended in .model.gz contain a model checkpoint file that stores information of all models tested
Files ended in .mldist contain site-specific subtitution rates determined by maximum likelihood (read with MS Excel or R with command: #tab=read.table('example.phy.mlrate',header=TRUE))
Files ended in .iqtree are the main report file with a textual representation of the final tree

#### From MAFFT alignment

GUC-Mm_N_Sanger2_aligned.fa.mldist
GUC-Mm_N_Sanger2_aligned.fa.treefile
GUC-Mm_N_Sanger2_aligned.fa.ckp.gz
GUC-Mm_N_Sanger2_aligned.fa.model.gz
GUC-Mm_N_Sanger2_aligned.fa.log
GUC-Mm_N_Sanger2_aligned.fa.iqtree
GUC-Mm_N_Sanger2_aligned.fa.contree
GUC-Mm_N_Sanger2_aligned.fa.bionj
GUC-Mm_N_Sanger2_aligned.fa.splits.nex

#### From MUSCLE alignment

GUC-Mm_N_Sanger2_MUSCLE.afa.iqtree
GUC-Mm_N_Sanger2_MUSCLE.afa.ckp.gz
GUC-Mm_N_Sanger2_MUSCLE.afa.treefile
GUC-Mm_N_Sanger2_MUSCLE.afa.contree
GUC-Mm_N_Sanger2_MUSCLE.afa.splits.nex
GUC-Mm_N_Sanger2_MUSCLE.afa.mldist
GUC-Mm_N_Sanger2_MUSCLE.afa.log
GUC-Mm_N_Sanger2_MUSCLE.afa.bionj

#### From MAFFT + GBLOCKS alignment

GUC-Mm_N_Sanger2_aligned.fa-gb.treefile
GUC-Mm_N_Sanger2_aligned.fa-gb.mldist
GUC-Mm_N_Sanger2_aligned.fa-gb.splits.nex
GUC-Mm_N_Sanger2_aligned.fa-gb.model.gz
GUC-Mm_N_Sanger2_aligned.fa-gb.log
GUC-Mm_N_Sanger2_aligned.fa-gb.contree
GUC-Mm_N_Sanger2_aligned.fa-gb.iqtree
GUC-Mm_N_Sanger2_aligned.fa-gb.ckp.gz
GUC-Mm_N_Sanger2_aligned.fa-gb.bionj

#### From MUSCLE + GBLOCKS alignment

GUC-Mm_N_Sanger2_MUSCLE.afa-gb.bionj
GUC-Mm_N_Sanger2_MUSCLE.afa-gb.iqtree
GUC-Mm_N_Sanger2_MUSCLE.afa-gb.log
GUC-Mm_N_Sanger2_MUSCLE.afa-gb.treefile
GUC-Mm_N_Sanger2_MUSCLE.afa-gb.ckp.gz
GUC-Mm_N_Sanger2_MUSCLE.afa-gb.contree
GUC-Mm_N_Sanger2_MUSCLE.afa-gb.splits.nex
GUC-Mm_N_Sanger2_MUSCLE.afa-gb.model.gz
GUC-Mm_N_Sanger2_MUSCLE.afa-gb.mldist

### RevBayes Bayesian topology inference from topology.Rev script

GUC-Mm_N_Sanger2_MUSCLE_18S.tre    Tree file in NEWICK format
GUC-Mm_N_Sanger2_MUSCLE_18S.log    Log file
GUC-Mm_N_Sanger2_MUSCLE_18S.trees    All 3001 marginal tree topologies 
GUC-Mm_N_Sanger2_MUSCLE_18S.cons.tree     Consensus tree file in NEWICK format with forced positive branch lengths


### RevBayes Time Tree analysis (scripts below)

GUC-Mm_N_Sanger2_MUSCLE_ML_rooted_multiD_ultra.tree    Constrain topology input file derived from GUC-Mm_N_Sanger2_MUSCLE.afa.contree pruned and ultrametrized in AccessionComparison_and_Backbone.R script
TimeTree_GUC-Mm_N_Sanger2_MUSCLE_mcmc.log    Log file
TimeTree_GUC-Mm_N_Sanger2_MUSCLE_mcmc.trees    All 3001 marginal tree branch length variants 
TimeTree_GUC-Mm_N_Sanger2_MUSCLE_mcmc_MAP.tre     MAP consensus tree file in NEWICK format (used in phyloMS_code.R for ancestral state analyses)

## Access other information

Links to other publicly accessible locations of the data:
  * https://github.com/antropoteuthis/salp_ecomorphology/tree/master/phylogeny/ML will contain the GTR+Gamma alternate versions of the ML inferences

Literature data on colony architecture was derived from: Madin, L. P. (1990). Aspects of jet propulsion in salps. Canadian Journal of Zoology, 68(4), 765-777.

Newly submitted GenBank accessions:

	Accession_number	Specimen
	OQ863569	28_Metcalfina_hexagona_D39-Mhex-B1       
	OQ863570	3_Ihlea_punctata_non-spotted_D24-Ipun-B-1
	OQ863571	4_Ihlea_punctata_D24-Ipun-B-1            
	OQ863572	13_Cyclosalpa_bakeri_D27-Cbak-B-1        
	OQ863573	6_Iasis_cylindrica_D22-Icyl-B-1          
	OQ863574	27_Cyclosalpa_polae_D38-Cpol-B-1         
	OQ863575	7_Ritteriella_amboinensis_D23-Ramb-B-1   
	OQ863576	26_Cyclosalpa_quadriluminis_D39-Cqua-OS-1
	OQ863577	31_Helicosalpa_virgula_D42-Hvir-B-1      
	OQ863578	9_Helicosalpa_virgula_MarcHughes_specimen2
	OQ863579	16_Ritteriella_amboinensis_D28-Ramb-OS-1 
	OQ863580	20_Ihlea_punctata_D31-Ipun-OS-1+2        
	OQ863581	12_Helicosalpa_younti_D27-Hyou-B-1       
	OQ863582	15_Cyclosalpa_pinnata_D28-Cpin-B-1       
	OQ863583	23_Ritteriella_retracta_D37-Rsp-B-1      
	OQ863584	24_Ritteriella_retracta_D37-Rret-OS-1     


## Code/Software

### R scripts (RStudio)

phyloMS_code.R    R code for the downstream analyses. Tree and data pruning, phylogenetic uncertainty assessment, SIMMAP and ML reconstructions of architectural characters, evolutionary model fitting and contrasts, phylogenetic signal estimation. 

Required packages:
	library(tidyverse)
	library(ape)
	library(phytools)
	library(reshape2)
	library(data.table)
	library(geiger)
	library(bayou)
	library(RColorBrewer)
	library(phylolm)

AccessionComparison_and_Backbone.R    R code for comparing the effects of including different GenBank accessions as ingroups and outgroups for the phylogenetic analysis; and to prune and ultrametrize the constrain ML tree for Byesian time tree analyses in RevBayes.

Required packages:
	require(ape)
	library(phytools)
	library(magrittr)

Alignment_ML.commands.sh.txt    List of all Shell commands used for phylogenetic inference including code to build alignments, GBLOCKS postprocessing, IQTree ML analyses, and Bayesian analyses.


### RevBayes scripts (run in RevBayes)

#### For Bayesian tree topology estimation

Authors: Sebastian Höhna, Michael Landis, Brian Moore, and Tracy Heath
Source: https://revbayes.github.io/tutorials/ctmc/

topology.Rev    Bayesian inference of phylogeny using a GTR+Gamma+Inv substitution model

#### For the relaxed molecular clock time tree estimation

Author: Tracy A. Heath
Source: https://revbayes.github.io/tutorials/clocks/
Input files: GUC-Mm+N+Sanger2_MUSCLE.afa (alignment), GUC-Mm_N_Sanger2_MUSCLE_ML_rooted_multiD_ultra.tree (tree file)

mcmc_TimeTree_salps.Rev    Time tree full model specification and MCMC set up for estimating time-calibrated phylogenies under a clock model.
m_BDP_salps.Rev    Time tree birth-death model specification file for fixed topology
m_UCLN_salps.Rev    Time tree UCLN relaxed-clock model specification file for branches
m_GTR_salps.Rev    Time tree General Time Reversible matrix model specification file

#### For categorical ancestral state reconstruction using ordered Markov models

Authors:  Sebastian Höhna, Will Freyman, April M. Wright, Michael J. Landis
Source: https://revbayes.github.io/tutorials/morph_ase/ase_free.html
Input files: morph.nex (tip states), GUCMmNSanger2_TimeTree_salp18Sphylo.tre (tree file)

mcmc_ase_freeK_RJ.Rev    Runs the full reversible-jump MCMC on an irreversible model with multiple rates for each developmental transition
mcmc_ase_singleRate.Rev    Runs the full reversible-jump MCMC on an irreversible model with a single rate parameter for all developmental transitions



