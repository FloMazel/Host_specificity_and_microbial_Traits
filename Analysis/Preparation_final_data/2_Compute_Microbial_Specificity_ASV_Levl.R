####################################################
####   Microbial specificity computations       ####
####################################################

rm(list=ls())
set.seed(1)
# Load packages
library(ape)
library(tidyverse)
library(vegan)
library(picante)
library(PhyloMeasures)
library(PerformanceAnalytics)

source("Scripts/Final/0_Utility_functions.R")

#######################################
####  Specificity at the ASV LEVEL ####
#######################################

null.model="richness" # conserve ASV prevalence and shuffle its distribution across host individuals

###################################
####  Song et al. dataset      ####
###################################

#metadata
SP_samples=read_csv("Data/Final_Processed_Song_Data/Subseted_Song_Metadata_withgroup_ID.csv")
table_Sp_count=table(SP_samples$Time_Tree_curated)

# read counts 
Bacteria_OTU_t_raref=read.table("Data/Final_Processed_Song_Data/ASVs_counts_filt_rar5000.txt",check.names = F)
dim(Bacteria_OTU_t_raref)
OTU_table=as.matrix(Bacteria_OTU_t_raref)
OTU_table[1:5,1:5]
OTU_tableN <- (OTU_table[,SP_samples$dada2_output_names])
colnames(OTU_tableN) <- as.character(SP_samples$SampleID) # change smaple names

# microbial taxonomy 
Bacterial_Taxo=read.table("Data/Final_Processed_Song_Data/ASVs_taxonomy_filt.txt",stringsAsFactors = F)
row.names(Bacterial_Taxo)=Bacterial_Taxo$ASV

# Host Phylogeny 
HostTree=read.tree("Data/EMP_data/total_timetree_names.bdtt-filt.nwk.tre")

# Subset samples so that each host species has at least 2 representative samples
# ------------------------------------------------------------------------------

subset_samples=subset(SP_samples,Time_Tree_curated%in%names(table_Sp_count)[table_Sp_count>1]) %>% 
  subset(Time_Tree_curated%in%HostTree$tip.label)

#Create host species * ASV table (filed with 0 and 1)
subset_samples$Pres=1
sample_to_host_species <- subset_samples %>% 
  select(SampleID,Time_Tree_curated,Pres) %>% 
  pivot_wider(names_from=Time_Tree_curated,values_from=Pres)

sample_to_host_species[is.na(sample_to_host_species)]=0
sample_to_host_speciesM <- as.matrix(sample_to_host_species[,-1])
row.names(sample_to_host_speciesM)=sample_to_host_species$SampleID

#Compute specificity
One_samples_5k=specificity(subset_name="SONG_morethan2_samples_per_host_5k",
                             Mammals_samples=subset_samples,
                             prevASV_table=OTU_tableN,
                             sample_to_host_species=sample_to_host_speciesM,
                             HostTree=HostTree,
                             Bacterial_Taxo=Bacterial_Taxo,
                              null.model=null.model)

# Subset samples so that each host species has 5 representative samples
# ------------------------------------------------------------------------------
subset_samples= SP_samples %>%
  subset(Time_Tree_curated%in%names(table_Sp_count)[table_Sp_count>5]) %>% 
  group_by(Time_Tree_curated) %>%
  sample_n(5,replace = F) %>%
  subset(Time_Tree_curated%in%HostTree$tip.label)

#Create host species * ASV table (filed with 0 and 1)
subset_samples$Pres=1
sample_to_host_species <- subset_samples %>% 
  select(SampleID,Time_Tree_curated,Pres) %>% 
  pivot_wider(names_from=Time_Tree_curated,values_from=Pres)

sample_to_host_species[is.na(sample_to_host_species)]=0
sample_to_host_speciesM <- as.matrix(sample_to_host_species[,-1])
row.names(sample_to_host_speciesM)=sample_to_host_species$SampleID

#Compute specificity
five_samples_5k=specificity(subset_name="SONG_5_samples_per_host_5k",
                             Mammals_samples=subset_samples,
                             prevASV_table=OTU_tableN,
                             sample_to_host_species=sample_to_host_speciesM,
                             HostTree=HostTree,
                             Bacterial_Taxo=Bacterial_Taxo,
                            null.model=null.model)


#Combine results and save
ASV_specificity=bind_cols(One_samples_5k,five_samples_5k[rownames(One_samples_5k),])
ASV_specificity  
dim(ASV_specificity)
write.table(ASV_specificity,"Data/Final_Processed_Song_Data/ASV_specificity_Song_5k.csv")


###############     ####################
####  Youngblut et al. dataset      ####
###############     ####################

# Load metadata (equivalence samples IDs - host ID)
Mammals_samples=read.table("Data/Youngblut_Data/Mammals_Updated_Metadata.txt",header=T,stringsAsFactors = F) %>% 
  mutate(TimeTree_returned=host_taxonomy_species,
         Time_Tree_curated=host_taxonomy_species,
         SampleID=sample_ID_2)

# microbial taxonomy 
Bacterial_Taxo=read.table("Data/Youngblut_Data/ASVs_taxonomy.txt",stringsAsFactors = F)

Bact_Tree=read.tree("Data/Youngblut_Data/FasTree_GTR_Gamma_mafft_aligned_ASVs_pfiltered_withouthUnclassified.tre")

rarOTUtableRAW=read.table("Data/Final_Processed_Youngblut_Data/ASVs_counts_filt_rar4900.txt")
# Rename the samples
rarOTUtableRAW = rarOTUtableRAW[Mammals_samples$ENA_Id,]
ENA=Mammals_samples$SampleID; names(ENA)=Mammals_samples$ENA_Id
rownames(rarOTUtableRAW)=ENA[rownames(rarOTUtableRAW)]

commonBact=intersect(colnames(rarOTUtableRAW),Bact_Tree$tip.label)

Bact_Tree=drop.tip(Bact_Tree,Bact_Tree$tip.label[!Bact_Tree$tip.label%in%commonBact])
rarOTUtable=t(rarOTUtableRAW[,commonBact])

table_Sp_count=table(Mammals_samples$Time_Tree_curated)
table_Sp_count

#Create host species * ASV table (filed with 0 and 1)
Mammals_samples$Value=1
sample_to_host_species=spread(Mammals_samples[,],key=TimeTree_returned,value=Value)

row.names(sample_to_host_species)=sample_to_host_species$SampleID
sample_to_host_species=sample_to_host_species[,colnames(sample_to_host_species)%in%Mammals_samples$TimeTree_returned]

sample_to_host_species[is.na(sample_to_host_species)]=0


# Host Phylogeny 
HostTree=read.tree("Data/Host_phylo_FaurbyExact.tre")
HostTree=drop.tip(HostTree,HostTree$tip.label[!HostTree$tip.label%in%Mammals_samples$host_taxonomy_species])
HostTree


##########################
#    Compute Specificty  #
##########################

#### 5K TABLE #####

# All host species with at least two samples 
subset_samples=subset(Mammals_samples,TimeTree_returned%in%names(table_Sp_count)[table_Sp_count>1]) # keep only samples in species with more than one sample
One_samples=specificity(subset_name="morethan2_samples_per_host_5k_Youngblut",
                        Mammals_samples=subset_samples,prevASV_table=rarOTUtable,
                        sample_to_host_species=sample_to_host_species,HostTree=HostTree,
                        Bacterial_Taxo=Bacterial_Taxo,
                        null.model=null.model)

#save
write.table(One_samples,"Data/Final_Processed_Youngblut_Data/ASV_specificity_5k_YoungBlut.csv")

