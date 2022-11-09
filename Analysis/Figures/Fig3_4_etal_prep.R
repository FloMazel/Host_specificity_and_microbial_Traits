### ------------------------------------------------------------ ###
##  Result : Variability in microbial specificity across microbes      ###
### -------------------------------------------------------------###

rm(list=ls())

#############################
######    Preparation  ######
#############################

## Load data and functions, define directory and ploting charatcteristics
##------------------------------------------------------------------------

# Packages 
library(ggplot2)
library(cowplot)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
source("Scripts/Submitted/0_Utility_functions.R")
library(cowplot)
library(ape)
library(phytools)


## colors for phylums
cols = c("Bacteroidetes"="#A6761D","Proteobacteria"="#E6AB02","Firmicutes"="#66A61E","Tenericutes"="#D95F02","Other" = "#666666")
Non_main_phyla <- c("Actinobacteria","Chlamydiae","Elusimicrobia","Epsilonbacteraeota", "Fibrobacteres","Fusobacteria", "Spirochaetes","Verrucomicrobia")


# Load data #
#############

# Specificity 
ASV_specificity_Song_raw=as_tibble(read.table("Data/Final_Processed_Song_Data/ASV_specificity_Song_5k.csv",stringsAsFactors = F))
ASV_specificityYoungblut_raw=as_tibble(read.table("Data/Final_Processed_Youngblut_Data/ASV_specificity_5k_YoungBlut.csv",stringsAsFactors = F))


Genus_Specificity=read.csv("Data/SubmissionData/Genus_Specificity_Song_Youngblut.csv",stringsAsFactors = F)


# FILTERING PARAMETERS #
########################

# Define thresholds for ASVs
N_sample_min_perASV = 2
N_reads_min_perASV = 100

# Define thresholds for Genera
N_ASV_min_perGenus = 0
q75 <- quantile(Genus_Specificity$observed_readCounts_SONG_morethan2_samples_per_host_5k_sum,.75,na.rm = T)
q75



# Subset DATA according to paremeters #
#######################################

ASV_specificity_Song_filt = ASV_specificity_Song_raw %>% 
  dplyr::rename(Genus=Genus_SONG_morethan2_samples_per_host_5k) %>% 
  dplyr::rename(Family=Family_SONG_morethan2_samples_per_host_5k)  %>% 
  dplyr::rename(Order=Order_SONG_morethan2_samples_per_host_5k)  %>% 
  dplyr::rename(Phylum=Phylum_SONG_morethan2_samples_per_host_5k)  %>% 
  subset(observed_samples_counts_SONG_morethan2_samples_per_host_5k>N_sample_min_perASV) %>% #keep only ASVs observed in at least 3 samples. 
  subset(observed_readCounts_SONG_morethan2_samples_per_host_5k>N_reads_min_perASV) %>% 
  subset(!is.infinite(sr.obs.z_SONG_morethan2_samples_per_host_5k)&!is.infinite(sr.obs.z_SONG_morethan2_samples_per_host_5k)) %>% 
  dplyr::mutate(Phylum_simplified= ifelse(Phylum %in% Non_main_phyla, "Other",Phylum))


# Yougblut
ASV_specificity_Y_filt =  ASV_specificityYoungblut_raw %>% 
  dplyr::rename(Genus=Genus_morethan2_samples_per_host_5k_Youngblut,
         Phylum = Phylum_morethan2_samples_per_host_5k_Youngblut) %>% 
  subset(observed_samples_counts_morethan2_samples_per_host_5k_Youngblut>N_sample_min_perASV) %>% #keep only ASVs observed in at least 3 samples. 
  subset(observed_readCounts_morethan2_samples_per_host_5k_Youngblut>N_reads_min_perASV) %>% 
  subset(!is.infinite(sr.obs.z_morethan2_samples_per_host_5k_Youngblut)&!is.infinite(sr.obs.z_morethan2_samples_per_host_5k_Youngblut)) %>% 
  mutate(Phylum_simplified= ifelse(Phylum %in% Non_main_phyla, "Other",Phylum))


# Genus level
Abgenius = Genus_Specificity %>% 
  subset(observed_readCounts_SONG_morethan2_samples_per_host_5k_sum > q75)  %>% 
  pull(Genus) %>% as.character()

Rich_Genus =  ASV_specificity_Song_filt %>% 
  group_by(Genus) %>% 
  summarise(n=n()) %>% 
  subset(n>N_ASV_min_perGenus) %>%
  na.omit() %>% 
  pull(Genus) %>% as.character()


Genus_Specificity = Genus_Specificity %>% 
  subset(Genus%in%Abgenius) %>% 
  subset(Genus%in%Rich_Genus)  %>% 
  mutate(Phylum_simplified = ifelse(Phylum %in% Non_main_phyla, "Other",Phylum)) %>% 
  mutate(Genus_Phylo = ifelse(Genus=="Hafnia-Obesumbacterium","Hafnia",
                              ifelse(Genus%in%c("Coprococcus_2","Coprococcus_3"),"Coprococcus",
                                     ifelse(Genus%in%c("Ruminiclostridium_5","Ruminiclostridium_6","Ruminiclostridium_9"),"Ruminiclostridium",
                                            ifelse(Genus%in%c("Ruminococcus_1","Ruminococcus_2"),"Ruminococcus",Genus)))))


ASV_specificity_Song_raw %>%
  mutate(Abundant=(observed_samples_counts_SONG_morethan2_samples_per_host_5k>N_sample_min_perASV)&(observed_readCounts_SONG_morethan2_samples_per_host_5k>N_reads_min_perASV)) %>% #keep only ASVs observed in at least 3 samples. 
  mutate(Assignation=Genus_SONG_morethan2_samples_per_host_5k%in%Abgenius) %>%
  group_by(Assignation,Abundant) %>% 
  summarise(nASVs=n(),
            nReads=sum(observed_readCounts_SONG_morethan2_samples_per_host_5k)) %>% 
  mutate(propASV=nASVs/sum(nASVs)/2,
         propReads=nReads/sum(nReads)/2)


ASV_specificity_Song <- ASV_specificity_Song_filt %>% # from 22000 ASvs (100% reads) to 4700 ASVs (92% reads) to 2000 ASVs (102 genera, for 50% reads)
  subset(Genus%in%Abgenius) %>% 
  subset(Genus%in%Rich_Genus) 


ASV_specificityYoungblut_raw %>%
  mutate(Abundant = (observed_samples_counts_morethan2_samples_per_host_5k_Youngblut>N_sample_min_perASV)&(observed_readCounts_morethan2_samples_per_host_5k_Youngblut>N_reads_min_perASV)) %>% #keep only ASVs observed in at least 3 samples. 
  mutate(Assignation=Genus_morethan2_samples_per_host_5k_Youngblut%in%Abgenius) %>%
  group_by(Assignation,Abundant) %>% 
  summarise(nASVs=n(),
            nReads=sum(observed_readCounts_morethan2_samples_per_host_5k_Youngblut),) %>% 
  mutate(propASV=nASVs/sum(nASVs)/2,
         propReads=nReads/sum(nReads)/2)


ASV_specificity_Y <- ASV_specificity_Y_filt %>% 
  subset(Genus%in%Abgenius) %>% 
  subset(Genus%in%Rich_Genus)


######################################
######################################
#    Microbial genus phylogeny       #
######################################
######################################

SSU_tree <- ape::read.tree ("Data/SSURefNR99_1200_slv_132_2.ntree")

# Rename tree tips to Acc ID
Tree_tips <- tibble(rawNames=SSU_tree$tip.label) %>% 
  separate(rawNames,sep=",",into=c("Acc"),remove = F) %>% 
  separate(Acc,sep="'",into=c(NA,"AccessionID"))

SSU_tree_Accname = SSU_tree
SSU_tree_Accname$tip.label = Tree_tips$AccessionID

# Select tips that correspond to our 102 genera and subset the tree
taxo_slv <- read_delim("Data/taxmap_slv_ssu_ref_nr_132.txt") %>% 
  separate(path,sep=";",into=as.character(1:6))
taxo_slv$Genus <- gsub(" ","_",taxo_slv$`6`)

taxo_slv_myGenus <- taxo_slv %>% 
  subset(Genus%in%Genus_Specificity$Genus)

SSU_tree_Accname_MyGenera <- keep.tip(SSU_tree_Accname,taxo_slv_myGenus$primaryAccession)
SSU_tree_Accname_MyGenera

myGenusTrees <- list()
ntree=50

for (i in 1:ntree)
{
  # Randomly subsample one tip per genus and produce a tree
  print(i)
  OneTip_list <- taxo_slv_myGenus %>% 
    group_by(Genus) %>% 
    sample_n(1) %>% 
    pull(primaryAccession)
  SSU_tree_Accname_MyGeneraOneTip <- keep.tip(SSU_tree_Accname_MyGenera,OneTip_list)
  
  # rename tree by genus name 
  #SSU_tree_Accname_MyGeneraOneTip$tip.label = taxo_slv_myGenus$Genus[match(OneTip_list,taxo_slv_myGenus$primaryAccession)]
  SSU_tree_Accname_MyGeneraOneTip$tip.label = taxo_slv_myGenus$Genus[match(SSU_tree_Accname_MyGeneraOneTip$tip.label,taxo_slv_myGenus$primaryAccession)]
  #plot(SSU_tree_Accname_MyGeneraOneTip)
  SSU_tree_Accname_MyGeneraOneTip$node.label<-NULL
  myGenusTrees[[i]] <-  SSU_tree_Accname_MyGeneraOneTip
}

saveRDS(myGenusTrees, "Data/MicrobialGenusTrees.RDS")

  
