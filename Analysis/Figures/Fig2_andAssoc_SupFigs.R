### ----------------------------------------------------------------------- ###
##  Result 1: Effect of host species on microbiome composition              ###
### ------------------------------------------------------------------------###

rm(list=ls())

## Load data and functions, define directory and ploting charatcteristics
##------------------------

# Load packages
library(reshape)
library(tidyverse)
library(parallel)
library(vegan)
library(ggConvexHull)
library(cowplot)
library(RColorBrewer)
source("Scripts/Submitted/0_Utility_functions.R")


#define list of data
Sample_list=list();Sample_list[["Song"]]=list();Sample_list[["Youngblut"]]=list()
Betas=list();Betas[["Song"]]=list();Betas[["Youngblut"]]=list()

# Load YOUNGBLUT metadata (equivalence samples IDs - host ID)
SP_samples=as_tibble(read.table("Data/Youngblut_Data/Mammals_Updated_Metadata.txt",header=T,stringsAsFactors = F)) %>%
  mutate(Time_Tree_curated=TimeTree_returned)

# Subset samples so that each host species has at least 2 representative samples
table_Sp_count=table(SP_samples$Time_Tree_curated)
Sample_list[["Youngblut"]][["Samples_from_species_with_more_than_1_sample"]] <- SP_samples %>% 
  subset(Time_Tree_curated%in%names(table_Sp_count)[table_Sp_count>1])

table_Sp_count=table(Sample_list[["Youngblut"]][["Samples_from_species_with_more_than_1_sample"]]$Time_Tree_curated)
mean(c(table_Sp_count))
sd(c(table_Sp_count))
summary(c(table_Sp_count))


# Load SONG metadata (equivalence samples IDs - host ID)
SP_samples=read_csv("Data/Final_Processed_Song_Data/Subseted_Song_Metadata_withgroup_ID.csv")
table_Sp_count=table(SP_samples$Time_Tree_curated)

# Subset samples so that each host species has at least 2 representative samples
Sample_list[["Song"]][["Samples_from_species_with_more_than_1_sample"]] <- SP_samples %>% 
  subset(Time_Tree_curated%in%names(table_Sp_count)[table_Sp_count>1])

table_Sp_count=table(Sample_list[["Song"]][["Samples_from_species_with_more_than_1_sample"]]$Time_Tree_curated)
mean(c(table_Sp_count))
sd(c(table_Sp_count))
summary(c(table_Sp_count))

# Keep species with at least 5 samples (subset the one with > 5 samples to 5 samples and keep all sampels for species with n=5 samples exactly) 
Sample_list[["Song"]][["Samples_from_species_with_5_samples"]] = SP_samples %>% 
  subset(Time_Tree_curated%in%names(table_Sp_count)[table_Sp_count>5]) %>% 
  group_by(Time_Tree_curated) %>%
  sample_n(5,replace = F)


# Number of host and samples
dim(Sample_list[["Song"]][["Samples_from_species_with_more_than_1_sample"]])
length(unique(Sample_list[["Song"]][["Samples_from_species_with_more_than_1_sample"]]$Time_Tree_curated))
dim(Sample_list[["Youngblut"]][["Samples_from_species_with_more_than_1_sample"]])
length(unique(Sample_list[["Youngblut"]][["Samples_from_species_with_more_than_1_sample"]]$Time_Tree_curated))

dim(Sample_list[["Song"]][["Samples_from_species_with_more_than_1_sample"]])+dim(Sample_list[["Youngblut"]][["Samples_from_species_with_more_than_1_sample"]])
length(unique(c(Sample_list[["Youngblut"]][["Samples_from_species_with_more_than_1_sample"]]$Time_Tree_curated,Sample_list[["Song"]][["Samples_from_species_with_more_than_1_sample"]]$Time_Tree_curated)))


# Load YOUNGBLUT Beta-diversity matrices 
BetaPath="Intermediate_Data/Beta_diversity/Youngblut_et_al/"
Betas[["Youngblut"]]=lapply(paste(BetaPath,list.files(BetaPath),sep=""),readRDS)
slices=sapply(list.files(BetaPath),FUN=function(x){strsplit(strsplit(x,split=".RDS")[[1]],"simi")[[1]][2]}) # Different phylogenetic resolution to defince microbial units. 
names(Betas[["Youngblut"]])=slices

# Load SONG Beta-diversity matrices 
BetaPath="Intermediate_Data/Beta_diversity/Song_et_al/"
Betas[["Song"]]=lapply(paste(BetaPath,list.files(BetaPath),sep=""),readRDS)
slices=sapply(list.files(BetaPath),FUN=function(x){strsplit(strsplit(x,split=".RDS")[[1]],"simi")[[1]][2]}) # Different phylogenetic resolution to defince microbial units. 
names(Betas[["Song"]])=slices

# Renames Song Beta matrix dada2
Betas[["Song"]][['0']]=Betas[["Song"]][['0']][,SP_samples$dada2_output_names,SP_samples$dada2_output_names]
dimnames(Betas[["Song"]][['0']])[2][[1]]=dimnames(Betas[["Song"]][['0']])[3][[1]]=as.character(SP_samples$SampleID)



##----------------------------##
##  Run the PERMANOVA models  ##
##----------------------------##

# Differenet models (srata-= diet, order or non)
sratas=c("Diet","Order","None")
#sratas=c("Order","None")

#Different beta-diversity metrics
metrics=c("Bray","Jac_TT")

# Two different datasets: samples from host species that have exactly 5 samples, or all samples with more than 1 sample. 
datasets_YoungBlut=c("Samples_from_species_with_more_than_1_sample")
datasets_Song=c("Samples_from_species_with_more_than_1_sample","Samples_from_species_with_5_samples")

# Create the combination of parameters 
Comb_Song=expand.grid(sratas,metrics,slices,datasets_Song,"Song",stringsAsFactors=F)
Comb_Youngblut=expand.grid(sratas,metrics,slices,datasets_YoungBlut,"Youngblut",stringsAsFactors=F)
Comb=rbind(Comb_Song,Comb_Youngblut)
Disp_results=tibble(Phylo_resolution=Comb$Var3,Sample_set=Comb$Var4,Metrics=Comb$Var2,Strata=Comb$Var1,Data_set=Comb$Var5,F_value=NA,P_value=NA,Number_of_Species=NA)

# Run the models
nb_param=dim(Disp_results)[1]

#Disp_results[2,]
wrapper_adonis1(1,Disp_results,Betas,Sample_list)
#Results=mclapply(1:3,try(wrapper_adonis1),Disp_results,Betas,Sample_list,mc.cores = 4)

Results=mclapply(1:nb_param,try(wrapper_adonis1),Disp_results,Betas,Sample_list,mc.cores = 4)
Results_staked=do.call(rbind,Results)
Disp_results$F_value=Results_staked[,1]
Disp_results$P_value=Results_staked[,2]
Disp_results$Number_of_Species=Results_staked[,4]
Disp_results$R2=Results_staked[,3]

#Disp_results <- read_csv("Outputs/Supp_Tables/Supp_Table1_Both_Datasets.csv")

write_csv(Disp_results,"Outputs/Supp_Tables/Supp_Table1_Both_Datasets_withgroup_ID.csv")

# Print result for the main text
Disp_results=read.csv("Outputs/Supp_Tables/Supp_Table1_Both_Datasets_withgroup_ID.csv",stringsAsFactors = F)
a = Disp_results %>%
  subset(Phylo_resolution=="0"&Metrics=="Bray"&Sample_set=="Samples_from_species_with_more_than_1_sample")

#Format the table for exporting (Supp. Table 1)
Disp_results <- Disp_results %>% 
  mutate(`Individuals choice` = changeSampleName(Sample_set),
          `Beta-diversity metric` = changeBetaName(Metrics),
          `Controling for (Permutation block)` = Strata,
         `Data set` = Data_set,
         `Number of species` = Number_of_Species,
         `pseudo-F value` = round(F_value,2),
         R2=round(R2,2),
         P_value=ifelse(P_value==.001,"<.001",P_value)
         ) %>% 
  select(`Data set`,`Number of species`,`Individuals choice`,`Beta-diversity metric`,`Controling for (Permutation block)`,
         `pseudo-F value`,R2,P_value)
         
write.csv(Disp_results,"Redaction/V2/Supp_data/Supp_Table_2.csv")

##------------------------------##
##  Run the PERMADISPERS models ##
##------------------------------##

# Differenet models (srata-= diet, order or non)
sratas=c("None")

#Different beta-diversity metrics
metrics=c("Bray" ,"Jac_TT")

# Two different datasets: samples from host species that have exactly 5 samples, or all samples with more than 1 sample. 
datasets_YoungBlut=c("Samples_from_species_with_more_than_1_sample")
datasets_Song=c("Samples_from_species_with_more_than_1_sample","Samples_from_species_with_5_samples")

# Create the combination of parameters 
Comb_Song=expand.grid(sratas,metrics,'0',datasets_Song,"Song",stringsAsFactors=F)
Comb_Youngblut=expand.grid(sratas,metrics,'0',datasets_YoungBlut,"Youngblut",stringsAsFactors=F)
Comb=rbind(Comb_Song,Comb_Youngblut)
Disp_results=tibble(Phylo_resolution=Comb$Var3,Sample_set=Comb$Var4,Metrics=Comb$Var2,Strata=Comb$Var1,Data_set=Comb$Var5,F_value=NA,P_value=NA,Number_of_Species=NA)

# Run the models
nb_param=dim(Disp_results)[1]

wrapper_betadisper(1,Disp_results,Betas,Sample_list)

Results=mclapply(1:nb_param,try(wrapper_betadisper),Disp_results,Betas,Sample_list,mc.cores = 4)

Results_staked=do.call(rbind,Results)
Disp_results$F_value=Results_staked[,1]
Disp_results$P_value=Results_staked[,2]
Disp_results$Number_of_Species=Results_staked[,3]

#Format the table for exporting (Supp. Table 1)
#write.csv(Disp_results,"Outputs/Tables/Supp_Table2_Both_DatasetsNEW.csv")
write.csv(Disp_results,"Redaction/V2/Supp_data/Supp_Table_3.csv")

# Print result for the main text
Disp_results=read.csv("Redaction/V2/Supp_data/Supp_Table_3.csv",stringsAsFactors = F)
Disp_results %>%
  subset(Phylo_resolution=="0"&Metrics=="Bray"&Sample_set=="Samples_from_species_with_more_than_1_sample")

Disp_results_clean <- Disp_results %>% 
  mutate(`Individuals choice` = changeSampleName(Sample_set),
         `Beta-diversity metric` = changeBetaName(Metrics),
         `Data set` = Data_set,
         `Number of species` = Number_of_Species,
         `pseudo-F value` = round(F_value,2),
         P_value=ifelse(P_value==.001,"<.001",P_value)
  ) %>% 
  select(`Data set`,`Number of species`,`Individuals choice`,`Beta-diversity metric`,
         `pseudo-F value`,P_value)

write.csv(Disp_results_clean,"Redaction/V2/Supp_data/Supp_Table_3.csv")


##------------##
##  Figure 1  ##
##------------##

# Subsample to 5 indiviudal for plots 
SP_samples=read_csv("Data/Final_Processed_Song_Data/Subseted_Song_Metadata_withgroup_ID.csv")

table_Sp_count=table(SP_samples$Time_Tree_curated)
subset_for_plots= SP_samples %>%
  subset(Time_Tree_curated%in%names(table_Sp_count)[table_Sp_count>4]) %>%
  group_by(Time_Tree_curated) %>% 
  sample_n(5)

length(unique(subset_for_plots$Time_Tree_curated))

## Extract beta-diversities values 
s=as.data.frame(subset_for_plots); row.names(s)=s$SampleID
Beta=Betas[["Song"]][["0"]]["Bray",s$SampleID,s$SampleID]

## NMDS plot
##-----------
set.seed(10)

# Run the NMDS
NMDS_run=metaMDS(Beta,k=2,try=30)
NMDS_run$stress

# Store the results
NMDS=data.frame(Axe1=NMDS_run$points[,1],Axe2=NMDS_run$points[,2],Species=s$Time_Tree_curated,Order=s$Order) %>%
  mutate(Taxonomy=ifelse(Order%in%c("Carnivora","Cetartiodactyla","Chiroptera","Perisodactyla","Primates","Rodentia"),as.character(Order),"Other"))

#Plot the results ggplot(NMDS_O,aes(y=Axe1,x=Axe2,group=Species,col=Order))+geom_point()+geom_convexhull(alpha=.01)+MyTheme
NMDS_plot= NMDS %>% 
  ggplot(aes(y=Axe1,x=Axe2,group=Species,col=Taxonomy,fill=Taxonomy)) +
  geom_point(size=.5)+geom_convexhull(alpha=0.05,size=.2)+
  xlab("NMDS axis 2")+ ylab("NMDS axis 1")+
  xlim(c(-.3,.5))+
  MyTheme + theme(legend.position = c(.9,.8))+  
  ggtitle("A") + 
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2")

NMDS_plot

## Beta diversity boxplot
##-----------------------

# reformat the beta-diversities 

Beta=Betas[["Song"]][["0"]]["Bray",s$SampleID,s$SampleID]
#Beta=Betas[["Song"]][["0"]]["Jac",,]
dimnames(Betas[["Song"]][["0"]])[[1]]
dim(Beta)

reshaped_Beta=as.data.frame(melt(Beta)) %>% 
  subset(!X1==X2)
reshaped_Beta$Sp_1=s[as.character(reshaped_Beta$X1),"Time_Tree_curated"]
reshaped_Beta$Sp_2=s[as.character(reshaped_Beta$X2),"Time_Tree_curated"]
reshaped_Beta=reshaped_Beta %>% mutate(comparison=ifelse(Sp_1==Sp_2,"Within","Between"))

head(reshaped_Beta)

taxo <- SP_samples %>% 
  dplyr::select(Time_Tree_curated,Taxonomy_Genus) %>% 
  group_by(Time_Tree_curated) %>% 
  summarise(Genus=unique(Taxonomy_Genus))

reshaped_Beta <- reshaped_Beta %>% 
  left_join(taxo,by=c("Sp_1"="Time_Tree_curated")) %>% 
  left_join(taxo,by=c("Sp_2"="Time_Tree_curated")) %>% 
  
  mutate(comparisonGenus=ifelse(Genus.x==Genus.y,"WithinGenus","BetweenGenus"))

head(reshaped_Beta)
  
  
summrayBeta = reshaped_Beta %>% 
  group_by(comparisonGenus) %>% 
  summarise(Mean = mean(value),
            Sd=sd(value),
            Median = median(value))

summrayBeta

summrayBeta = reshaped_Beta %>% 
  group_by(comparison) %>% 
  summarise(Mean = mean(value),
            Median = median(value))

(1-0.995)*100

Boxplot_beta <- reshaped_Beta %>% 
  ggplot(aes(y=1-value,x=comparison))+
  MyTheme+geom_boxplot(alpha=.01,outlier.size=.1)+
  ylab("1 - Bray Curtis dissimilarity")+xlab("Within and Between Host species comparisons")+
  ggtitle("B")
Boxplot_beta



# Final plot + export
fig_1=plot_grid(NMDS_plot,Boxplot_beta, align = "h", nrow = 1, rel_widths = c(2,1))
fig_1

# Figure size (mm)
units_fig_size="mm"
typical_h =89
two_column_w = 183
one_column_w = 89 

ggsave(filename = "Submission/Mol_Ecol/V2/Main_figures/Fig_2.pdf",
       plot = fig_1,
       device = "pdf",
       height = typical_h,
       width = two_column_w,
       units=units_fig_size)

##-----------------##
##  SUpp Figure 1 (Youngblut)  ##
##-----------------##

# Subsample to 5 indiviudal for plots 
SP_samples=as_tibble(read.table("Data/Youngblut_Data/Mammals_Updated_Metadata.txt",header=T,stringsAsFactors = F)) %>%
  mutate(Time_Tree_curated=TimeTree_returned)

table_Sp_count=table(SP_samples$Time_Tree_curated)
subset_for_plots= SP_samples %>%
  subset(Time_Tree_curated%in%names(table_Sp_count)[table_Sp_count>2])# %>%
  #group_by(Time_Tree_curated) %>% 
  #sample_n(5)

subset_for_plots %>% 
  summarise(n=n(),
            nsp=n_distinct(Time_Tree_curated))

## Extract beta-diversities values 
s=as.data.frame(subset_for_plots); row.names(s)=s$SampleID
Beta=Betas[["Youngblut"]][["0"]]["Bray",s$SampleID,s$SampleID]

## NMDS plot
##-----------

# Run the NMDS
NMDS_run=metaMDS(Beta,k=2)

# Store the results
NMDS=data.frame(Axe1=NMDS_run$points[,1],Axe2=NMDS_run$points[,2],Species=s$Time_Tree_curated,Order=s$Order) %>%
  mutate(Taxonomy=ifelse(Order%in%c("Carnivora","Cetartiodactyla","Chiroptera","Perisodactyla","Primates","Rodentia"),as.character(Order),"Other"))

#Plot the results ggplot(NMDS_O,aes(y=Axe1,x=Axe2,group=Species,col=Order))+geom_point()+geom_convexhull(alpha=.01)+MyTheme
NMDS_plot=ggplot(NMDS,aes(y=Axe1,x=Axe2,group=Species,col=Taxonomy))+geom_point()+geom_convexhull(alpha=.05)+
  MyTheme+  ggtitle("A") + scale_color_brewer(palette="Dark2")
NMDS_plot

## Beta diversity boxplot
##-----------------------

# reformat the beta-diversities 
reshaped_Beta=as.data.frame(melt(Beta)) %>% 
  subset(!X1==X2)

reshaped_Beta$Sp_1=s[as.character(reshaped_Beta$X1),"Time_Tree_curated"]
reshaped_Beta$Sp_2=s[as.character(reshaped_Beta$X2),"Time_Tree_curated"]
reshaped_Beta=reshaped_Beta %>% mutate(comparison=ifelse(Sp_1==Sp_2,"Within","Between"))

Boxplot_beta=ggplot(reshaped_Beta,aes(y=1-value,x=comparison))+MyTheme+geom_boxplot()+
  ylab("1 - Bray Curtis dissimilarity")+xlab("Within and Between Host species comparisons")+
  ggtitle("B")


# Final plot + export ! 
fig_1=plot_grid(NMDS_plot,Boxplot_beta, align = "h", nrow = 1, rel_widths = c(2,1))
fig_1
ggsave(filename = "Redaction/V2/Supp_data/Supp_Figure_1.pdf",
       plot = fig_1,
       device = "pdf",
       height = typical_h,
       width = two_column_w,
       units=units_fig_size)




