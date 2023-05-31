#### Packages ----
library(ape)
library(phylotools)
library(doBy)
library(tidyverse)
library(phangorn)

library(MCMCglmm)
library(metafor)
library(dplyr)
library(reshape2)
library(phytools)

## Import tree into R  ----
tree_newick <- read.tree("Original_tree.txt") # In newick format

write.nexus(tree_newick, file="Original_tree.nex")

tree_nexus <- read.nexus("Original_tree.nex")
summary(tree_nexus)

tree <- tree_nexus

is.ultrametric(tree)
tree <- chronoMPL(tree) # Make it ultrametric if it isn't

species_data <- read.csv("pangenome_species.csv", header=T)

## Make data frame of tip labels of tree
#write.csv(data.frame(treeSp = tree$tip.label), file = "treeSp.csv")

#write.csv(data.frame(treeSp = dataTree$tip.label), file = "treeSp_dataTree.csv")

#### Replacing with proper tip name
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Rhodospirillales_Acetobacteraceae_Acetobacter_pasteurianus_386B")] <- "Acetobacter_pasteurianus"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Pseudomonadales_Moraxellaceae_Acinetobacter_baumannii_AB30")] <- "Acinetobacter_baumannii"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Aeromonadales_Aeromonadaceae_Aeromonas_hydrophila_AL09_71")] <- "Aeromonas_hydrophila"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Alteromonadales_Alteromonadaceae_Alteromonas_macleodii_Black_Sea_11")] <- "Alteromonas_mediterranea" #same species, different name
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Bacillales_Bacillaceae_Bacillus_anthracis_52_G")] <- "Bacillus_anthracis"
tree$tip.label[which(tree$tip.label == "Bacteria_Actinobacteria_Actinobacteria_Actinobacteridae_Bifidobacteriales_Bifidobacteriaceae_Bifidobacterium_animalis_animalis_ATCC_25527")] <- "Bifidobacterium_animalis" 
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Rhizobiales_Brucellaceae_Brucella_abortus_A13334")] <- "Brucella_abortus" 
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Buchnera_Buchnera_aphidicola_str._APS_Acyrthosiphon_pisum")] <- "Buchnera_aphidicola"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Betaproteobacteria_Burkholderiales_Burkholderiaceae_Burkholderia_cenocepacia_AU_1054")] <- "Burkholderia_cenocepacia"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_deltaepsilon_subdivisions_Epsilonproteobacteria_Campylobacterales_Campylobacteraceae_Campylobacter_jejuni_jejuni_NCTC_11168")] <- "Campylobacter_jejuni"	
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiaceae_Clostridium_difficile_CF5")] <- "Clostridioides_difficile"	#same species
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Legionellales_Coxiellaceae_Coxiella_burnetii_CbuK_Q154")] <- "Coxiella_burnetii"	
tree$tip.label[which(tree$tip.label == "Bacteria_Actinobacteria_Actinobacteria_Actinobacteridae_Actinomycetales_Propionibacterineae_Propionibacteriaceae_Propionibacterium_acnes_KPA171202")] <- "Cutibacterium_acnes"	#same species
tree$tip.label[which(tree$tip.label == "Bacteria_Chloroflexi_Dehalococcoidetes_Dehalococcoidales_Dehalococcoidaceae_Dehalococcoides_mccartyi_CBDB1")] <- "Dehalococcoides_mccartyi"	
tree$tip.label[which(tree$tip.label == "Bacteria_Bacteroidetes_Chlorobi_group_Bacteroidetes_Flavobacteriia_Flavobacteriales_Flavobacteriaceae_Elizabethkingia_anophelis_NUH6")] <- "Elizabethkingia_anophelis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Enterobacter_cloacae_SCF1")] <- "Enterobacter_cloacae"	
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Lactobacillales_Enterococcaceae_Enterococcus_faecalis_D32")] <- "Enterococcus_faecalis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Escherichia_coli_K12__W3110")] <- "Escherichia_coli"	
tree$tip.label[which(tree$tip.label == "Bacteria_Fusobacteria_Fusobacteriia_Fusobacteriales_Fusobacteriaceae_Fusobacterium_nucleatum_nucleatum_ATCC_25586")] <- "Fusobacterium_nucleatum"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Rhodospirillales_Acetobacteraceae_Granulibacter_bethesdensis_CGDNIH1")] <- "Granulibacter_bethesdensis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Pasteurellales_Pasteurellaceae_Haemophilus_influenzae_KR494")] <- "Haemophilus_influenzae"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Klebsiella_oxytoca_HKOPL1")] <- "Klebsiella_oxytoca"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Pasteurellales_Pasteurellaceae_Mannheimia_haemolytica_D174")] <- "Mannheimia_haemolytica"	
tree$tip.label[which(tree$tip.label == "Bacteria_Actinobacteria_Actinobacteria_Actinobacteridae_Actinomycetales_Corynebacterineae_Mycobacteriaceae_Mycobacterium_abscessus_bolletii_50594")] <- "Mycobacteroides_abscessus"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Betaproteobacteria_Neisseriales_Neisseriaceae_Neisseria_gonorrhoeae_FA_1090")] <- "Neisseria_gonorrhoeae"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Pasteurellales_Pasteurellaceae_Pasteurella_multocida_43137")] <- "Pasteurella_multocida"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Thiotrichales_Piscirickettsiaceae_Piscirickettsia_salmonis_LF_89")] <- "Piscirickettsia_salmonis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Bacteroidetes_Chlorobi_group_Bacteroidetes_Bacteroidia_Bacteroidales_Porphyromonadaceae_Porphyromonas_gingivalis_W83")] <- "Porphyromonas_gingivalis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Cyanobacteria_Prochlorales_Prochlorococcaceae_Prochlorococcus_marinus_str._MIT_9301")] <- "Prochlorococcus_marinus"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Proteus_mirabilis_WGLW4")] <- "Proteus_mirabilis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Pseudomonadales_Pseudomonadaceae_Pseudomonas_aeruginosa_B136_33")] <- "Pseudomonas_aeruginosa"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Betaproteobacteria_Burkholderiales_Burkholderiaceae_Ralstonia_solanacearum_GMI1000")] <- "Ralstonia_solanacearum"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Rhizobiales_Rhizobiaceae_RhizobiumAgrobacterium_group_Rhizobium_leguminosarum_bv._trifolii_CB782")] <- "Rhizobium_leguminosarum"	
tree$tip.label[which(tree$tip.label == "Bacteria_Bacteroidetes_Chlorobi_group_Bacteroidetes_Flavobacteriia_Flavobacteriales_Flavobacteriaceae_Riemerella_anatipestifer_RA_CH_2")] <- "Riemerella_anatipestifer"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Salmonella_enterica_Enterica_sv._Heidelberg_41578")] <- "Salmonella_enterica"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Rhizobiales_Rhizobiaceae_SinorhizobiumEnsifer_group_Ensifer_meliloti_GR4")] <- "Sinorhizobium_meliloti"	
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Bacillales_Staphylococcaceae_Staphylococcus_aureus_502A")] <- "Staphylococcus_aureus"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Xanthomonadales_Xanthomonadaceae_Stenotrophomonas_maltophilia_K279a")] <- "Stenotrophomonas_maltophilia"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Vibrionales_Vibrionaceae_Vibrio_anguillarum_775_I")] <- "Vibrio_anguillarum" #checked genus phylogeny
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Xanthomonadales_Xanthomonadaceae_Xanthomonas_oryzae_ATCC_35933")] <- "Xanthomonas_oryzae" #checked genus phylogeny
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Xanthomonadales_Xanthomonadaceae_Xylella_fastidiosa_sandyi_Ann_1")] <- "Xylella_fastidiosa"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Yersinia_pestis_Nepal516")] <- "Yersinia_pestis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Sphingomonadales_Sphingomonadaceae_Zymomonas_mobilis_mobilis_CP4")] <- "Zymomonas_mobilis"	


#### Add species to tree, replace other species in same genus
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Betaproteobacteria_Burkholderiales_Alcaligenaceae_Bordetella_bronchiseptica_RB50")] <- "Bordetella_holmesii"	
tree$tip.label[which(tree$tip.label == "Bacteria_Chlamydiae_Verrucomicrobia_group_Chlamydiae_Chlamydiia_Chlamydiales_Chlamydiaceae_ChlamydiaChlamydophila_group_Chlamydophila_abortus_S263")] <- "Chlamydia_psittaci" #checked genus phylogeny
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Citrobacter_rodentium_ICC168")] <- "Citrobacter_freundii"
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiaceae_Clostridium_acetobutylicum_EA_2018")] <- "Clostridium_botulinum" #checked genus phylogeny
tree$tip.label[which(tree$tip.label == "Bacteria_Actinobacteria_Actinobacteria_Actinobacteridae_Actinomycetales_Corynebacterineae_Corynebacteriaceae_Corynebacterium_argentoratense_DSM_44202")] <- "Corynebacterium_glutamicum" #checked genus phylogeny
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Thiotrichales_Francisellaceae_Francisella_philomiragia_subsp._philomiragia_ATCC_25017")] <- "Francisella_tularensis" #checked genus phylogeny
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_deltaepsilon_subdivisions_Epsilonproteobacteria_Campylobacterales_Helicobacteraceae_Helicobacter_cetorum_MIT_99_5656")] <- "Helicobacter_pylori"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Legionellales_Legionellaceae_Legionella_micdadei_ATCC_33218")] <- "Legionella_pneumophila"
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Lactobacillales_Leuconostocaceae_Leuconostoc_kimchii_IMSNU11154")] <- "Leuconostoc_mesenteroides"
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Bacillales_Listeriaceae_Listeria_innocua_sv._6a_Clip11262")] <- "Listeria_monocytogenes"
tree$tip.label[which(tree$tip.label == "Bacteria_Tenericutes_Mollicutes_Mycoplasmatales_Mycoplasmataceae_Mycoplasma_agalactiae")] <- "Mycoplasma_bovis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Bacillales_Paenibacillaceae_Paenibacillus_lactis_154")] <- "Paenibacillus_polymyxa"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Rhodobacterales_Rhodobacteraceae_Phaeobacter_gallaeciensis_2.10")] <- "Phaeobacter_inhibens"
tree$tip.label[which(tree$tip.label == "Bacteria_Bacteroidetes_Chlorobi_group_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotella_dentalis_ES_2772_DSM_3688")] <- "Prevotella_intermedia"
tree$tip.label[which(tree$tip.label == "Bacteria_Actinobacteria_Actinobacteria_Propionibacteriales_Propionibacteriaceae_Propionibacterium_CG1_02_FULL_Propionibacterium_60_36")] <- "Propionibacterium_freudenreichii"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Rickettsiales_Rickettsiaceae_Rickettsieae_Rickettsia_conorii_Malish_7")] <- "Rickettsia_prowazekii"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Serratia_fonticola_RB_25")] <- "Serratia_marcescens"
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Lactobacillales_Streptococcaceae_Streptococcus_mutans_GS_5")] <- "Streptococcus_thermophilus"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Xanthomonadales_Xanthomonadaceae_Xanthomonas_axonopodis_Xac29_1")] <- "Xanthomonas_campestris" #checked genus phylogeny
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Lactobacillales_Streptococcaceae_Lactococcus_garvieae_Lg2")] <- "Lactococcus_lactis"
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Lactobacillales_Lactobacillaceae_Lactobacillus_acidophilus_NCFM")] <- "Lactobacillus_delbrueckii" #in most basal clade, random between delbrueckii and helveticus which was chosen

#### REMOVE EVERYTHING ELSE FROM TREE

# Drop tips that aren't in data set
# Creates what you don't want
dropTip2<-tree$tip.label[which(is.na(match(tree$tip.label, species_data$species)))]
# Gets rid of it
dataTree<-drop.tip(tree, dropTip2, trim.internal=T)
plot(dataTree)

# Check which species aren't in tree
species_data$species[which((species_data$species %in% dataTree$tip.label)==FALSE)]


#### ADDING SPECIES TO GENUS
## A. pittii
which(dataTree$tip.label=="Acinetobacter_baumannii")		# 17 (edge number of this tip)
which(dataTree$edge[,2] == 17)	  #  43 (row in edge matrix with this end node)
dataTree$edge.length[43]									# 0.32 (edge length)
dataTree <- bind.tip(dataTree, "Acinetobacter_pittii", edge.length=NULL, where=17, position=dataTree$edge.length[43]/4) #4 if adding in first species, /2 if clustering further

## Bordetella pertussis
which(dataTree$tip.label=="Bordetella_holmesii")		# 27 (edge number of this tip)
which(dataTree$edge[,2] == 27)	  #  61 (row in edge matrix with this end node)
dataTree$edge.length[61]									# 0.36 (edge length)
dataTree <- bind.tip(dataTree, "Bordetella_pertussis", edge.length=NULL, where=27, position=dataTree$edge.length[61]/4) #4 if adding in first species, /2 if clustering further

## E. hormachei
which(dataTree$tip.label=="Enterobacter_cloacae")		# 4 (edge number of this tip)
which(dataTree$edge[,2] == 4)	  #  25 (row in edge matrix with this end node)
dataTree$edge.length[25]									# 0.022 (edge length)
dataTree <- bind.tip(dataTree, "Enterobacter_hormaechei", edge.length=NULL, where=4, position=dataTree$edge.length[25]/4) 

## E. faecium
which(dataTree$tip.label=="Enterococcus_faecalis")		# 57 (edge number of this tip)
which(dataTree$edge[,2] == 57)	  #  120 (row in edge matrix with this end node)
dataTree$edge.length[120]									# 0.15(edge length)
dataTree <- bind.tip(dataTree, "Enterococcus_faecium", edge.length=NULL, where=57, position=dataTree$edge.length[120]/4)

## Fusobacterium periodonticum
which(dataTree$tip.label=="Fusobacterium_nucleatum")		# 48 (edge number of this tip)
which(dataTree$edge[,2] == 48)	  #  96 (row in edge matrix with this end node)
dataTree$edge.length[96]									# 1.08(edge length)
dataTree <- bind.tip(dataTree, "Fusobacterium_periodonticum", edge.length=NULL, where=48, position=dataTree$edge.length[96]/8) #very long branch, do /8 to make more similar

## Neisseria meningitidis
which(dataTree$tip.label=="Neisseria_gonorrhoeae")		# 32 (edge number of this tip)
which(dataTree$edge[,2] == 32)	  #  69 (row in edge matrix with this end node)
dataTree$edge.length[69]									# 0.4 (edge length)
dataTree <- bind.tip(dataTree, "Neisseria_meningitidis", edge.length=NULL, where=32, position=dataTree$edge.length[69]/4) 

##R. phaseoli
which(dataTree$tip.label=="Rhizobium_leguminosarum")		# 34 (edge number of this tip)
which(dataTree$edge[,2] == 34)	  #  78 (row in edge matrix with this end node)
dataTree$edge.length[78]									# 0.064 (edge length)
dataTree <- bind.tip(dataTree, "Rhizobium_phaseoli", edge.length=NULL, where=34, position=dataTree$edge.length[78]/4) 

##S.epidermidis
which(dataTree$tip.label=="Staphylococcus_aureus")		# 65 (edge number of this tip)
which(dataTree$edge[,2] == 65)	  #  133 (row in edge matrix with this end node)
dataTree$edge.length[133]									# 0.375 (edge length)
dataTree <- bind.tip(dataTree, "Staphylococcus_epidermidis", edge.length=NULL, where=65, position=dataTree$edge.length[133]/4)


# Check which species aren't in tree
species_data$species[which((species_data$species %in% dataTree$tip.label)==FALSE)]

plot(dataTree)
#### Adding multiple species to  genus

## BACILLUS
##B. pumilus
which(dataTree$tip.label=="Bacillus_anthracis")		# 67 (edge number of this tip)
which(dataTree$edge[,2] == 67)	  #  136 (row in edge matrix with this end node)
dataTree$edge.length[136]									# 0.375 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_pumilus", edge.length=NULL, where=67, position=dataTree$edge.length[136]/4)

##B. subtilis
which(dataTree$tip.label=="Bacillus_anthracis")		# 67 (edge number of this tip)
which(dataTree$edge[,2] == 67)	  #  137 (row in edge matrix with this end node)
dataTree$edge.length[137]									# 0.09 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_subtilis", edge.length=NULL, where=67, position=dataTree$edge.length[137]/2) # /2 when clustering

##B. cereus
which(dataTree$tip.label=="Bacillus_anthracis")		# 67 (edge number of this tip)
which(dataTree$edge[,2] == 67)	  #  138 (row in edge matrix with this end node)
dataTree$edge.length[138]									# 0.04 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_cereus", edge.length=NULL, where=67, position=dataTree$edge.length[138]/2) 

##B. thuringiensis
which(dataTree$tip.label=="Bacillus_cereus")		# 68 (edge number of this tip)
which(dataTree$edge[,2] == 68)	  #  140 (row in edge matrix with this end node)
dataTree$edge.length[140]									# 0.02 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_thuringiensis", edge.length=NULL, where=68, position=dataTree$edge.length[140]/2) 

##B. licheniformis
which(dataTree$tip.label=="Bacillus_subtilis")		# 70 (edge number of this tip)
which(dataTree$edge[,2] == 70)	  #  143 (row in edge matrix with this end node)
dataTree$edge.length[143]									# 0.04 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_licheniformis", edge.length=NULL, where=70, position=dataTree$edge.length[143]/2) 

##B. velezensis
which(dataTree$tip.label=="Bacillus_licheniformis")		# 71 (edge number of this tip)
which(dataTree$edge[,2] == 71)	  #  145 (row in edge matrix with this end node)
dataTree$edge.length[145]									# 0.02 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_velezensis", edge.length=NULL, where=71, position=dataTree$edge.length[145]/2) 

##B. amyloliquefaciens
which(dataTree$tip.label=="Bacillus_velezensis")		# 72 (edge number of this tip)
which(dataTree$edge[,2] == 72)	  #  147 (row in edge matrix with this end node)
dataTree$edge.length[147]									# 0.01 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_amyloliquefaciens", edge.length=NULL, where=72, position=dataTree$edge.length[147]/2) 

plot(dataTree)


## BIFIDOBACTERIUM
##B. breve
which(dataTree$tip.label=="Bifidobacterium_animalis")		# 52 (edge number of this tip)
which(dataTree$edge[,2] == 52)	  #  107 (row in edge matrix with this end node)
dataTree$edge.length[107]									# 0.385 (edge length)
dataTree <- bind.tip(dataTree, "Bifidobacterium_breve", edge.length=NULL, where=52, position=dataTree$edge.length[107]/4) 

##B. longum
which(dataTree$tip.label=="Bifidobacterium_breve")		# 53 (edge number of this tip)
which(dataTree$edge[,2] == 53)	  #  109 (row in edge matrix with this end node)
dataTree$edge.length[109]									# 0.096 (edge length)
dataTree <- bind.tip(dataTree, "Bifidobacterium_longum", edge.length=NULL, where=53, position=dataTree$edge.length[109]/2) 

plot(dataTree)

## BRUCELLA
##B. suis
which(dataTree$tip.label=="Brucella_abortus")		# 37 (edge number of this tip)
which(dataTree$edge[,2] == 37)	  #  82 (row in edge matrix with this end node)
dataTree$edge.length[82]									# 0.19 (edge length)
dataTree <- bind.tip(dataTree, "Brucella_suis", edge.length=NULL, where=37, position=dataTree$edge.length[82]/4) 

##B. melitensis
which(dataTree$tip.label=="Brucella_abortus")		# 37 (edge number of this tip)
which(dataTree$edge[,2] == 37)	  #  83 (row in edge matrix with this end node)
dataTree$edge.length[83]									# 0.04 (edge length)
dataTree <- bind.tip(dataTree, "Brucella_melitensis", edge.length=NULL, where=37, position=dataTree$edge.length[83]/2) 

plot(dataTree)

## BURKHOLDERIA
##B. thailandensis
which(dataTree$tip.label=="Burkholderia_cenocepacia")		# 31 (edge number of this tip)
which(dataTree$edge[,2] == 31)	  #  68 (row in edge matrix with this end node)
dataTree$edge.length[68]									# 0.15 (edge length)
dataTree <- bind.tip(dataTree, "Burkholderia_thailandensis", edge.length=NULL, where=31, position=dataTree$edge.length[68]/4) 

##B. mallei
which(dataTree$tip.label=="Burkholderia_thailandensis")		# 32 (edge number of this tip)
which(dataTree$edge[,2] == 32)	  #  70 (row in edge matrix with this end node)
dataTree$edge.length[70]									# 0.03 (edge length)
dataTree <- bind.tip(dataTree, "Burkholderia_mallei", edge.length=NULL, where=32, position=dataTree$edge.length[70]/2) 

##B. pseudomallei
which(dataTree$tip.label=="Burkholderia_mallei")		# 33 (edge number of this tip)
which(dataTree$edge[,2] == 33)	  #  72 (row in edge matrix with this end node)
dataTree$edge.length[72]									# 0.01 (edge length)
dataTree <- bind.tip(dataTree, "Burkholderia_pseudomallei", edge.length=NULL, where=33, position=dataTree$edge.length[72]/2) 

plot(dataTree)

## CAMPYLOBACTER
##C. fetus
which(dataTree$tip.label=="Campylobacter_jejuni")		# 49 (edge number of this tip)
which(dataTree$edge[,2] == 49)	  #  101 (row in edge matrix with this end node)
dataTree$edge.length[101]									# 0.27 (edge length)
dataTree <- bind.tip(dataTree, "Campylobacter_fetus", edge.length=NULL, where=49, position=dataTree$edge.length[101]/4) 

##C. coli
which(dataTree$tip.label=="Campylobacter_jejuni")		# 49 (edge number of this tip)
which(dataTree$edge[,2] == 49)	  #  102 (row in edge matrix with this end node)
dataTree$edge.length[102]									# 0.06 (edge length)
dataTree <- bind.tip(dataTree, "Campylobacter_coli", edge.length=NULL, where=49, position=dataTree$edge.length[102]/2) 

plot(dataTree)

## CHLAMYDIA
##C. trachomatis
which(dataTree$tip.label=="Chlamydia_psittaci")		# 52 (edge number of this tip)
which(dataTree$edge[,2] == 52)	  #  106 (row in edge matrix with this end node)
dataTree$edge.length[106]									# 1.05 (edge length)
dataTree <- bind.tip(dataTree, "Chlamydia_trachomatis", edge.length=NULL, where=52, position=dataTree$edge.length[106]/8) ##very long branch, so /8

##C. muridarum
which(dataTree$tip.label=="Chlamydia_trachomatis")		# 53 (edge number of this tip)
which(dataTree$edge[,2] == 53)	  #  108 (row in edge matrix with this end node)
dataTree$edge.length[108]									# 0.13 (edge length)
dataTree <- bind.tip(dataTree, "Chlamydia_muridarum", edge.length=NULL, where=53, position=dataTree$edge.length[108]/2) 

plot(dataTree)

## CORYENBACTERIUM
##C. diptheriae
which(dataTree$tip.label=="Corynebacterium_glutamicum")		# 67 (edge number of this tip)
which(dataTree$edge[,2] == 67)	  #  135 (row in edge matrix with this end node)
dataTree$edge.length[135]									# 0.18 (edge length)
dataTree <- bind.tip(dataTree, "Corynebacterium_diphtheriae", edge.length=NULL, where=67, position=dataTree$edge.length[135]/4) 

##C. pseudotuberculosis
which(dataTree$tip.label=="Corynebacterium_diphtheriae")		# 68 (edge number of this tip)
which(dataTree$edge[,2] == 68)	  #  137 (row in edge matrix with this end node)
dataTree$edge.length[137]									# 0.04 (edge length)
dataTree <- bind.tip(dataTree, "Corynebacterium_pseudotuberculosis", edge.length=NULL, where=68, position=dataTree$edge.length[137]/2) 

##C. ulcerans
which(dataTree$tip.label=="Corynebacterium_pseudotuberculosis")		# 69 (edge number of this tip)
which(dataTree$edge[,2] == 69)	  #  139 (row in edge matrix with this end node)
dataTree$edge.length[139]									# 0.02 (edge length)
dataTree <- bind.tip(dataTree, "Corynebacterium_ulcerans", edge.length=NULL, where=69, position=dataTree$edge.length[139]/2) 

plot(dataTree)

## KLEBSIELLA
## K. aerogenes
which(dataTree$tip.label=="Klebsiella_oxytoca")		# 6 (edge number of this tip)
which(dataTree$edge[,2] == 6)	  #  28 (row in edge matrix with this end node)
dataTree$edge.length[28]									# 0.026 (edge length)
dataTree <- bind.tip(dataTree, "Klebsiella_aerogenes", edge.length=NULL, where=6, position=dataTree$edge.length[28]/4) 

## Klebsiella_pneumoniae
which(dataTree$tip.label=="Klebsiella_aerogenes")		# 7 (edge number of this tip)
which(dataTree$edge[,2] == 7)	  #  30 (row in edge matrix with this end node)
dataTree$edge.length[30]									# 0.007 (edge length)
dataTree <- bind.tip(dataTree, "Klebsiella_pneumoniae", edge.length=NULL, where=7, position=dataTree$edge.length[30]/2) # divide by 2 as within genus should be equal branching 

plot(dataTree)

### LACTOBACILLUS
##L.sakei
which(dataTree$tip.label=="Lactobacillus_delbrueckii")		# 75 (edge number of this tip)
which(dataTree$edge[,2] == 75)	  #  156 (row in edge matrix with this end node)
dataTree$edge.length[156]									# 0.31 (edge length)
dataTree <- bind.tip(dataTree, "Lactobacillus_sakei", edge.length=NULL, where=75, position=dataTree$edge.length[156]/4) 

##L. paracasei
which(dataTree$tip.label=="Lactobacillus_sakei")		# 76 (edge number of this tip)
which(dataTree$edge[,2] == 76)	  #  158 (row in edge matrix with this end node)
dataTree$edge.length[158]									# 0.07 (edge length)
dataTree <- bind.tip(dataTree, "Lactobacillus_paracasei", edge.length=NULL, where=76, position=dataTree$edge.length[158]/2) 

##L. casei
which(dataTree$tip.label=="Lactobacillus_paracasei")		# 77 (edge number of this tip)
which(dataTree$edge[,2] == 77)	  #  160 (row in edge matrix with this end node)
dataTree$edge.length[160]									# 0.03 (edge length)
dataTree <- bind.tip(dataTree, "Lactobacillus_casei", edge.length=NULL, where=77, position=dataTree$edge.length[160]/2) 

##L. rhamnosus
which(dataTree$tip.label=="Lactobacillus_casei")		# 78 (edge number of this tip)
which(dataTree$edge[,2] == 78)	  #  160 (row in edge matrix with this end node)
dataTree$edge.length[162]									# 0.019 (edge length)
dataTree <- bind.tip(dataTree, "Lactobacillus_rhamnosus", edge.length=NULL, where=78, position=dataTree$edge.length[162]/2) 

##L. plantarum
which(dataTree$tip.label=="Lactobacillus_sakei")		# 76 (edge number of this tip)
which(dataTree$edge[,2] == 76)	  #  159 (row in edge matrix with this end node)
dataTree$edge.length[159]									# 0.03 (edge length)
dataTree <- bind.tip(dataTree, "Lactobacillus_plantarum", edge.length=NULL, where=76, position=dataTree$edge.length[159]/2) 

##L. fermentum
which(dataTree$tip.label=="Lactobacillus_plantarum")		# 77 (edge number of this tip)
which(dataTree$edge[,2] == 77)	  #  161 (row in edge matrix with this end node)
dataTree$edge.length[161]									# 0.019 (edge length)
dataTree <- bind.tip(dataTree, "Lactobacillus_fermentum", edge.length=NULL, where=77, position=dataTree$edge.length[161]/2) 

##L. brevis
which(dataTree$tip.label=="Lactobacillus_fermentum")		# 78 (edge number of this tip)
which(dataTree$edge[,2] == 78)	  #  163 (row in edge matrix with this end node)
dataTree$edge.length[163]									# 0.009 (edge length)
dataTree <- bind.tip(dataTree, "Lactobacillus_brevis", edge.length=NULL, where=78, position=dataTree$edge.length[163]/2) 

##L. helveticus
which(dataTree$tip.label=="Lactobacillus_delbrueckii")		# 75 (edge number of this tip)
which(dataTree$edge[,2] == 75)	  #  157 (row in edge matrix with this end node)
dataTree$edge.length[157]									# 0.07 (edge length)
dataTree <- bind.tip(dataTree, "Lactobacillus_helveticus", edge.length=NULL, where=75, position=dataTree$edge.length[157]/2) 

plot(dataTree)

## MYCOBACTERIUM (INC MYCOBACTEROIDES)
##M. avium
which(dataTree$tip.label=="Mycobacteroides_abscessus")		# 68 (edge number of this tip)
which(dataTree$edge[,2] == 68)	  #  138 (row in edge matrix with this end node)
dataTree$edge.length[138]									# 0.18 (edge length)
dataTree <- bind.tip(dataTree, "Mycobacterium_avium", edge.length=NULL, where=68, position=dataTree$edge.length[138]/4) 

##M. avium
which(dataTree$tip.label=="Mycobacterium_avium")		# 69 (edge number of this tip)
which(dataTree$edge[,2] == 69)	  #  140 (row in edge matrix with this end node)
dataTree$edge.length[140]									# 0.04 (edge length)
dataTree <- bind.tip(dataTree, "Mycobacterium_tuberculosis", edge.length=NULL, where=69, position=dataTree$edge.length[140]/2) 

##M. bovis
which(dataTree$tip.label=="Mycobacterium_tuberculosis")		# 70 (edge number of this tip)
which(dataTree$edge[,2] == 70)	  #  142 (row in edge matrix with this end node)
dataTree$edge.length[142]									# 0.02 (edge length)
dataTree <- bind.tip(dataTree, "Mycobacterium_bovis", edge.length=NULL, where=70, position=dataTree$edge.length[142]/2) 

plot(dataTree)

## MYCOPLASMA
##M. gallisepticum
which(dataTree$tip.label=="Mycoplasma_bovis")		# 102 (edge number of this tip)
which(dataTree$edge[,2] == 102)	  #  205 (row in edge matrix with this end node)
dataTree$edge.length[205]									# 0.49 (edge length)
dataTree <- bind.tip(dataTree, "Mycoplasma_gallisepticum", edge.length=NULL, where=102, position=dataTree$edge.length[205]/4) 

##M. mycoides
which(dataTree$tip.label=="Mycoplasma_bovis")		# 102 (edge number of this tip)
which(dataTree$edge[,2] == 102)	  #  206 (row in edge matrix with this end node)
dataTree$edge.length[206]									# 0.122 (edge length)
dataTree <- bind.tip(dataTree, "Mycoplasma_mycoides", edge.length=NULL, where=102, position=dataTree$edge.length[206]/2) 

##M. pneumoniae
which(dataTree$tip.label=="Mycoplasma_gallisepticum")		# 104 (edge number of this tip)
which(dataTree$edge[,2] == 104)	  #  209 (row in edge matrix with this end node)
dataTree$edge.length[209]									# 0.122 (edge length)
dataTree <- bind.tip(dataTree, "Mycoplasma_pneumoniae", edge.length=NULL, where=104, position=dataTree$edge.length[209]/2) 

plot(dataTree)

## PSEUDOMONAS
## P. putida
which(dataTree$tip.label=="Pseudomonas_aeruginosa")		# 19 (edge number of this tip)
which(dataTree$edge[,2] == 19)	  #  48 (row in edge matrix with this end node)
dataTree$edge.length[48]									# 0.32 (edge length)
dataTree <- bind.tip(dataTree, "Pseudomonas_putida", edge.length=NULL, where=19, position=dataTree$edge.length[48]/4) 

## P. fluorescens
which(dataTree$tip.label=="Pseudomonas_putida")		# 20 (edge number of this tip)
which(dataTree$edge[,2] == 20)	  #  50 (row in edge matrix with this end node)
dataTree$edge.length[50]									# 0.08 (edge length)
dataTree <- bind.tip(dataTree, "Pseudomonas_fluorescens", edge.length=NULL, where=20, position=dataTree$edge.length[50]/2) 

## P. syringae
which(dataTree$tip.label=="Pseudomonas_fluorescens")		# 21 (edge number of this tip)
which(dataTree$edge[,2] == 21)	  #  52 (row in edge matrix with this end node)
dataTree$edge.length[52]									# 0.04 (edge length)
dataTree <- bind.tip(dataTree, "Pseudomonas_syringae", edge.length=NULL, where=21, position=dataTree$edge.length[52]/2) 

## P. stutzeri
which(dataTree$tip.label=="Pseudomonas_aeruginosa")		# 19 (edge number of this tip)
which(dataTree$edge[,2] == 19)	  #  49 (row in edge matrix with this end node)
dataTree$edge.length[49]									# 0.08 (edge length)
dataTree <- bind.tip(dataTree, "Pseudomonas_stutzeri", edge.length=NULL, where=19, position=dataTree$edge.length[49]/2) 

plot(dataTree)

## RICKETTSIA
## R. rickettsii
which(dataTree$tip.label=="Rickettsia_prowazekii")		# 53 (edge number of this tip)
which(dataTree$edge[,2] == 53)	  #  110 (row in edge matrix with this end node)
dataTree$edge.length[110]									# 0.72 (edge length)
dataTree <- bind.tip(dataTree, "Rickettsia_rickettsii", edge.length=NULL, where=53, position=dataTree$edge.length[110]/6) #Long so /6, but not as long as the ones with /8

## R. japonica
which(dataTree$tip.label=="Rickettsia_rickettsii")		# 54 (edge number of this tip)
which(dataTree$edge[,2] == 54)	  #  112 (row in edge matrix with this end node)
dataTree$edge.length[112]									# 0.12 (edge length)
dataTree <- bind.tip(dataTree, "Rickettsia_japonica", edge.length=NULL, where=54, position=dataTree$edge.length[112]/2) 

plot(dataTree)

## STREPTOCOCCUS
## S. suis
which(dataTree$tip.label=="Streptococcus_thermophilus")		# 96 (edge number of this tip)
which(dataTree$edge[,2] == 96)	  #  197 (row in edge matrix with this end node)
dataTree$edge.length[197]									# 0.117 (edge length)
dataTree <- bind.tip(dataTree, "Streptococcus_suis", edge.length=NULL, where=96, position=dataTree$edge.length[197]/4)

## S. pyogenes
which(dataTree$tip.label=="Streptococcus_thermophilus")		# 96 (edge number of this tip)
which(dataTree$edge[,2] == 96)	  #  198 (row in edge matrix with this end node)
dataTree$edge.length[198]									# 0.029 (edge length)
dataTree <- bind.tip(dataTree, "Streptococcus_pyogenes", edge.length=NULL, where=96, position=dataTree$edge.length[198]/2)

## S. agalactiae
which(dataTree$tip.label=="Streptococcus_pyogenes")		# 97 (edge number of this tip)
which(dataTree$edge[,2] == 97)	  #  200 (row in edge matrix with this end node)
dataTree$edge.length[200]									# 0.014 (edge length)
dataTree <- bind.tip(dataTree, "Streptococcus_agalactiae", edge.length=NULL, where=97, position=dataTree$edge.length[200]/2)

## S. agalactiae
which(dataTree$tip.label=="Streptococcus_suis")		# 99 (edge number of this tip)
which(dataTree$edge[,2] == 99)	  #  203 (row in edge matrix with this end node)
dataTree$edge.length[203]									# 0.029 (edge length)
dataTree <- bind.tip(dataTree, "Streptococcus_pneumoniae", edge.length=NULL, where=99, position=dataTree$edge.length[203]/2)

plot(dataTree)

## VIBRIO
## V. vulnificus
which(dataTree$tip.label=="Vibrio_anguillarum")		# 16 (edge number of this tip)
which(dataTree$edge[,2] == 16)	  #  44 (row in edge matrix with this end node)
dataTree$edge.length[44]									# 0.21 (edge length)
dataTree <- bind.tip(dataTree, "Vibrio_vulnificus", edge.length=NULL, where=16, position=dataTree$edge.length[44]/4)

## V. alginolyticus
which(dataTree$tip.label=="Vibrio_vulnificus")		# 17 (edge number of this tip)
which(dataTree$edge[,2] == 17)	  #  46 (row in edge matrix with this end node)
dataTree$edge.length[46]									# 0.05 (edge length)
dataTree <- bind.tip(dataTree, "Vibrio_alginolyticus", edge.length=NULL, where=17, position=dataTree$edge.length[46]/2)

## V. cholerae
which(dataTree$tip.label=="Vibrio_vulnificus")		# 17 (edge number of this tip)
which(dataTree$edge[,2] == 17)	  #  47 (row in edge matrix with this end node)
dataTree$edge.length[47]									# 0.027 (edge length)
dataTree <- bind.tip(dataTree, "Vibrio_cholerae", edge.length=NULL, where=17, position=dataTree$edge.length[47]/2)

## V. parahaemolyticus
which(dataTree$tip.label=="Vibrio_alginolyticus")		# 19 (edge number of this tip)
which(dataTree$edge[,2] == 19)	  #  50 (row in edge matrix with this end node)
dataTree$edge.length[50]									# 0.027 (edge length)
dataTree <- bind.tip(dataTree, "Vibrio_parahaemolyticus", edge.length=NULL, where=19, position=dataTree$edge.length[50]/2)

plot(dataTree)

## XANTHOMONAS
## X. citri
which(dataTree$tip.label=="Xanthomonas_oryzae")		# 35 (edge number of this tip)
which(dataTree$edge[,2] == 35)	  #  78 (row in edge matrix with this end node)
dataTree$edge.length[78]									# 0.01 (edge length)
dataTree <- bind.tip(dataTree, "Xanthomonas_citri", edge.length=NULL, where=35, position=dataTree$edge.length[78]/2)

plot(dataTree)

## YERSINIA
## Yersinia_enterocolitica 
which(dataTree$tip.label=="Yersinia_pestis")		# 10 (edge number of this tip)
which(dataTree$edge[,2] == 10)	  #  36 (row in edge matrix with this end node)
dataTree$edge.length[36]									# 0.07 (edge length)
dataTree <- bind.tip(dataTree, "Yersinia_enterocolitica", edge.length=NULL, where=10, position=dataTree$edge.length[36]/4) 

## Yersinia_pseudotuberculosis
which(dataTree$tip.label=="Yersinia_pestis")		# 10 (edge number of this tip)
which(dataTree$edge[,2] == 10)	  #  37 (row in edge matrix with this end node)
dataTree$edge.length[37]									# 0.018 (edge length)
dataTree <- bind.tip(dataTree, "Yersinia_pseudotuberculosis", edge.length=NULL, where=10, position=dataTree$edge.length[37]/2) 

plot(dataTree)

# Check which species aren't in tree
species_data$species[which((species_data$species %in% dataTree$tip.label)==FALSE)]

## Remove 0 branch lengths by changing to 1 then remaking the tree ultrametric - repeat until no 0 branch lengths
dataTree$edge.length[dataTree$edge.length<=0]<-1
dataTree<-chronoMPL(dataTree)
which(dataTree$edge.length<=0)

summary(dataTree)
is.rooted(dataTree)

write.nexus(dataTree, file="panX_tree.nex")


