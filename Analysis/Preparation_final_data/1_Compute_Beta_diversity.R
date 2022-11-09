##-----------------------------##
## Beta-diversity computations ##
##-----------------------------##

rm(list=ls())

source("~/Documents/GitHub/BDTT/BDTT_functions.R")

library(vegan)
library(betapart)
library(abind)
library(Matrix)
library(ape)
library(tidyverse)


########## SONG DATA ################
# -------------------------------- #

# Load data 
#Bacteria_OTU_t_raref=read.table("/Users/fmazel/Data/Mammalian_microbiome_EMP/Dada2_output_separated_studies/ASVs_counts_rarefied5000.txt")
#dim(Bacteria_OTU_t_raref)
Bacteria_OTU_t_raref=read.table("Data/Final_Processed_Song_Data/ASVs_counts_filt_rar5000.txt",check.names = F)

dim(Bacteria_OTU_t_raref)
OTU_table=t(as.matrix(Bacteria_OTU_t_raref))
dim(OTU_table)
OTU_table[1:5,1:5]
#SP_samples=read_csv("Data/Song_Data/Metadata/Subseted_Song_Metadata.csv")
SP_samples=read_csv("Data/Song_Data/Metadata/Subseted_Song_Metadata_withgroup_ID.csv")

# Run BDTT 
#slices=rev(seq(from=0.1,to=1,by=0.1))
slices=0
for (i in slices)
{
  print(i)
  test=getBDTT(similarity=i,tree=NA,sampleOTUs=OTU_table,onlyBeta=T)
  saveRDS(test,file=paste("Data/Final_Processed_Song_Data/Beta/Song_BDTT_values5k_simi",i,".RDS",sep=""))
}


# On a different OTU table
Bacteria_OTU_t_raref=read.table("Data/Final_Processed_Song_Data/ASVs_counts_filt_rar5000_NoRare.txt",check.names = F)
dim(Bacteria_OTU_t_raref)
OTU_table=t(as.matrix(Bacteria_OTU_t_raref))
test=getBDTT(similarity=0,tree=NA,sampleOTUs=OTU_table,onlyBeta=T)
saveRDS(test,file=paste("Data/Final_Processed_Song_Data/Beta_NoRare/Song_BDTT_values5k_simi",i,".RDS",sep=""))



#"Intermediate_Data/Beta_diversity/Song_et_al/Song_BDTT_values5k_simi"

########## YOUNGBLUT et al.  DATA ################
# ---------------------------------------------- #

# Load data 
Bact_Tree=read.tree("Data/Youngblut_Data/FasTree_GTR_Gamma_mafft_aligned_ASVs_pfiltered_withouthUnclassified.tre")
#Asv_taxonomy=read.table("Data/Pre_processed/Dada2_ASVs/ASVs_taxonomy.txt",stringsAsFactors = F); Asv_taxonomy$ASV=rownames(Asv_taxonomy)
Mammals_samples=read.table("Data/Youngblut_Data/Mammals_Updated_Metadata.txt",header=T,stringsAsFactors = F)


#rarOTUtableRAW=read.table("Data/Youngblut_Data/rar4999_ASVs_counts.txt")
rarOTUtableRAW=read.table("Data/Final_Processed_Youngblut_Data/ASVs_counts_filt_rar4900.txt")

# Rename the samples
rarOTUtableRAW = rarOTUtableRAW[Mammals_samples$ENA_Id,]
ENA=Mammals_samples$SampleID; names(ENA)=Mammals_samples$ENA_Id
rownames(rarOTUtableRAW)=ENA[rownames(rarOTUtableRAW)]

Sample_distrib=table(table(MetadataMammals$host_taxonomy_species))
sum(Sample_distrib[!names(Sample_distrib)==1])
dim(rarOTUtableRAW)
Bact_Tree

commonBact=intersect(colnames(rarOTUtableRAW),Bact_Tree$tip.label)

Bact_Tree=drop.tip(Bact_Tree,Bact_Tree$tip.label[!Bact_Tree$tip.label%in%commonBact])
rarOTUtable=(rarOTUtableRAW[,commonBact])

slices=rev(seq(from=0,to=1,by=0.1))
slices=0
for (i in slices)
{
  print(i)
  test=getBDTT(similarity=i,tree=Bact_Tree,sampleOTUs=rarOTUtable,onlyBeta=T)
  saveRDS(test,file=paste("Data/Final_Processed_Youngblut_Data/BDTT_values_YoungBlut4900_simi",i,".RDS",sep=""))
}

#"Prep_Data/Beta_diversity_matrices/Youngblut_et_al/BDTT_values_YoungBlut4999_simi"
