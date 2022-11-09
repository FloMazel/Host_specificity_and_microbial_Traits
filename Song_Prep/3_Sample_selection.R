
######################################################
#####   Select samples from Song et al dataset  ######
######################################################

#rm(list=ls())
library(tidyverse)
library(purrr)
source("Scripts/Final/0_Utility_functions.R")

## ---------------------------------------------------------------------------------- ##
## Load samples with DNA and samples with METADATA and join the corresponding tables 
## ---------------------------------------------------------------------------------- ##

# Load Metadata
Metadata_online=read_csv("Data/Song_Data/Metadata/inline-supplementary-material-1.csv")
Song_samples=Metadata_online %>%
  subset(Taxonomy_Class=="Mammalia") %>% 
  mutate(TimeTree_returned=ifelse(Species_name=="Hipposideros camerunensis","Hipposideros_camerunensis",TimeTree_returned)) # small correction 



# Load sample list for which we have DNA data
Bacteria_OTU_t_raref=t(read.table("Data/Final_Processed_Song_Data/ASVs_counts_filt_rar5000.txt",check.names=F))
sample_names=rownames(Bacteria_OTU_t_raref)
Names=map_dfr(sample_names,re_shape_name)

# Get Basics info of those DNA data (e.g. number of reads, length of reads)
#raw_fastq_dir="/Users/fmazel/Data/Mammalian_microbiome_EMP/Raw_fastq/"; list_studies=list.files(raw_fastq_dir)
#All_samples <- makeSampleList(raw_fastq_dir) #make a list of all samples available in the dowloaded data
#fastq_info<-map2_dfr(All_samples$studyID,All_samples$sampleID,import_analyze_fastq,raw_fastq_dir) #extract read information (length, counts)
#write_csv(fastq_info,"Data/Final_Processed_Song_Data/Raw_fastq_files_Reads_characteristics.csv")

fastq_info=read_csv("Data/Final_Processed_Song_Data/Raw_fastq_files_Reads_characteristics.csv") %>%
  mutate(initial_fastq_name=sample_name) %>%
  select(-sample_name)

# Link DNA samples names with metadata names (new variable created: sampleID is the link with the metadata file)
Joined_names=full_join(fastq_info,Names,by='initial_fastq_name') %>%
  #mutate(duplicated_sample=ifelse(grepl('_duplicated',SampleID),T,F)) %>% # note the duplicated Song samples
  mutate(SampleID=str_split(initial_fastq_name,'.R1.gz',simplify = T)[,1]) %>%
  mutate(SampleID=str_split(SampleID,'.gz',simplify = T)[,1]) %>%
  mutate(SampleID=str_split(SampleID,'.R1_duplicated',simplify = T)[,1]) %>%
  mutate(SampleID=str_split(SampleID,'.fastq',simplify = T)[,1]) %>%
  mutate(SampleID=ifelse(studyID=="10285",paste("10285",SampleID,sep="."),SampleID)) %>% # rename Sanders samples (from MG_rast) 
  mutate(SampleID=ifelse(studyID=="1056",str_split(initial_fastq_name,'.653',simplify = T)[,1],SampleID)) %>%# rename Delsuc samples 
  mutate(SampleID=ifelse(studyID=="1056",gsub("seqs_","1056.",SampleID),SampleID)) %>%  # rename Delsuc samples 
  mutate(duplicated_sample=duplicated(SampleID) | duplicated(SampleID, fromLast=TRUE)) # note the duplicated Song samples


## ----------------------- ##
## SOME STATS on DNA DATA
## ----------------------- ##

Joined_names %>%
  subset(duplicated_sample==T) %>%
  pull(SampleID) %>%
  table() %>%
  table() # 162 samples are duplicated (150 or 250pb)

# Some overview
Joined_names %>%
  ggplot(aes(y=medianLength,x=duplicated_sample)) + geom_boxplot() +geom_jitter()

table(table(Joined_names$SampleID))

dupli = Joined_names %>% 
  subset(duplicated(SampleID) | duplicated(SampleID, fromLast=TRUE))# %>% 

## ------------------ ##
## SAMPLE SELECTION  
## ------------------ ##

# ----------------------------------------------------------------------------------------
# Steps of selection/modification

# 1. Remove duplicated samples names (keep the sample with highest number of reads) - could be because of multipole sequencing 
# 2. Select only samples that come form different individuals 
# 2bis. Remove deseased and new born/old individuals  
# 3. Host taxonomic cleaning 
# 4. Add Diet Metadata to the final table 
# ----------------------------------------------------------------------------------------

# -----------------------------
# 1. Duplicated sample names 
# -----------------------------

# Remove duplicata samplesID by choose one of the sample with highest number of reads 
Selected_Samples= Joined_names %>% 
  group_by(SampleID) %>% 
  top_n(1, read_number)

# --------------------------------------------------------
# 2. Select only samples that come form different individuals 
# --------------------------------------------------------

# 2.1. Duplicated individuals or groups of individuals 
# ----------------------------------------------------

# 2.1.1 Individuals information  
# -----------------------------

sum(is.na(Song_samples$individual_ID)) # there are lot of sampels without info on individuals (189)

# Stats 
Individual_ID_count = Song_samples %>% 
  subset(!is.na(individual_ID)) %>%
  group_by(individual_ID) %>%
  summarise(n_samples_per_indv=n())
table(Individual_ID_count$n_samples_per_indv)

# Duplicated and non_duplicated individual names 
duplicated_indv = Individual_ID_count  %>% 
  subset(n_samples_per_indv>1) %>%
  pull(individual_ID)

non_duplicated_indv = Individual_ID_count  %>% 
  subset(n_samples_per_indv==1) %>%
  pull(individual_ID)

length(duplicated_indv) # 387 indv with more than one sample. 

Duplicated_individual_data= Song_samples %>% 
  subset(individual_ID%in%duplicated_indv)

Summary_Individual_samples <- Duplicated_individual_data %>% 
  group_by(individual_ID) %>% 
  summarise(preservative=paste(unique(preservative),collapse = "_"),
            n_preservative=n_distinct(preservative),
            studyID=paste(unique(studyID),collapse="_"),
            country=paste(unique(country),collapse="_"),
            Order=Taxonomy_Order[1],
            TimeTree_returned=TimeTree_returned[1],
            Sample_size=n())

table(Summary_Individual_samples$n_preservative) # only one preservative by individual

# 2.1.2. Individuals groups information  
# -----------------------------

# stats on samples with a group ID entry 
group_ID_count = Song_samples %>% 
  subset(is.na(individual_ID)) %>%
  subset(!is.na(group_ID)) %>% 
  group_by(group_ID) %>%
  summarise(n_samples_per_gp=n())

table(group_ID_count$n_samples_per_gp)

# Does samples with individual ID also have a group ID? 
group_and_indv_ID = Song_samples %>% 
  subset(!is.na(individual_ID)) %>%
  subset(!is.na(group_ID)) 
# only 6 samples

duplicated_group_ID = group_ID_count  %>% 
  subset(n_samples_per_gp>1) %>%
  pull(group_ID)

non_duplicated_group_ID = group_ID_count  %>% 
  subset(n_samples_per_gp==1) %>%
  pull(group_ID)

Duplicated_group_ID= Song_samples %>% 
  subset(group_ID%in%duplicated_group_ID)

Summary_Group_samples<- Duplicated_group_ID %>% 
  group_by(group_ID) %>% 
  summarise(preservative=paste(unique(preservative),collapse = "_"),
            n_preservative=n_distinct(preservative),
            studyID=paste(unique(studyID),collapse="_"),
            country=paste(unique(country),collapse="_"),
            Order=Taxonomy_Order[1],
            TimeTree_returned=TimeTree_returned[1],
            Sample_size=n())

table(Summary_Group_samples$n_preservative) # only one preservative by individual


# 2.2. Sample selection 
# ---------------------

# We keep one sample per individual when individual information is available. 
# When a sample has no info on individual, but information on a grou of individual, we select samples coming
# from diferent group of individuals (in this case we are sure they also come from different individuals) 

Duplicated_Indv_choice  <- Duplicated_individual_data %>%
  group_by(individual_ID) %>% 
  sample_n(1) %>%  #choose randomly ONE sample per individual 
  pull(SampleID)

Duplicated_group_ID  <- Duplicated_group_ID %>%
  group_by(group_ID) %>% 
  sample_n(1) %>%  #choose randomly ONE sample per group of individuals 
  pull(SampleID)

# Final selection 
Final_metadata <- Song_samples %>% 
  subset(individual_ID%in%non_duplicated_indv|
           SampleID%in%Duplicated_Indv_choice|
           group_ID%in%non_duplicated_group_ID|
           SampleID%in%Duplicated_group_ID)

a = Final_metadata  %>% 
  subset(Taxonomy_Order=="Diprotodontia")

# Join with sample that have DNA 
Joined_Metadata <- Final_metadata %>% 
  inner_join(select(Selected_Samples,-studyID),by="SampleID") %>%
  subset(!is.na(dada2_output_names))  #only keep sampels that have been filtered (other have been excluded because of low number of reads)

# --------------------------------------------------------
# 2BIS.  Remove some individuals with too diverging life history
# --------------------------------------------------------

table(Joined_Metadata$LIFE_STAGE)
table(Joined_Metadata$healthy)

LifeStage_to_remove = c("newborn","old","juvenile")
health_to_remove = c("diseased","injured","roadkilled")

Joined_Metadata <- Joined_Metadata %>% 
  subset(!LIFE_STAGE%in%LifeStage_to_remove) %>% 
  subset(!healthy%in%health_to_remove) %>% 
  subset(!is.na(country)) %>% 
  subset(!is.na(preservative))

# Harmonize names 
Joined_Metadata <- Joined_Metadata %>%   
  mutate(country=ifelse(country=="United States of America","USA",country)) %>% # harmonize names 
  mutate(preservative=ifelse(preservative=="fta","FTA",preservative)) %>% # harmonize names 
  mutate(LIFE_STAGE=ifelse(LIFE_STAGE=="","unknown",LIFE_STAGE)) %>%
  mutate(sex=ifelse(sex%in%c("","F?","M?","Unknown","unknown"),NA,sex)) %>% 
  mutate(sex=ifelse(sex%in%c("F]"),"F",sex)) %>% 
  mutate(sex=ifelse(sex%in%c("Male_","F?"),unknown,sex)) %>% 
  mutate(sample_type=ifelse(sample_type=="feces","fecal",sample_type))


apply(Joined_Metadata,2,function(x){sum(is.na(x))})


# --------------------------------
# 3. Host taxonomic cleaning 
# --------------------------------

table(Joined_Metadata$TimeTree_returned)

# 3.1. Merge sub-species into one species 
# --------------------------------------

#Create a list of sub-species in the dataset
Time_Tree_subspecies=list(
  c("Miniopterus_natalensis_arenarius","Miniopterus_natalensis"),
  c("Elephas_maximus_indicus","Elephas_maximus"),
  c("Elephas_maximus_maximus","Elephas_maximus"),
  c("Equus_burchellii_chapmani","Equus_burchellii"),
  c("Equus_burchellii_quagga","Equus_burchellii"),
  c("Giraffa_camelopardalis_reticulata","Giraffa_camelopardalis"),
  c("Diceros_bicornis_michaeli","Diceros_bicornis"),
  c("Colobus_guereza_kikuyuensis","Colobus_guereza"),
  c("Gorilla_gorilla_gorilla","Gorilla_gorilla"),
  c("Pongo_pygmaeus_pygmaeus","Pongo_pygmaeus"),
  c("Saimiri_boliviensis_boliviensis","Saimiri_boliviensis"))

Time_Tree_subspecies <- data.frame(do.call(rbind,Time_Tree_subspecies))
colnames(Time_Tree_subspecies)=c("Sub_species_names","Species_name")

#Create a new column with updated species names 
Joined_Metadata <- Joined_Metadata %>% 
  mutate(is_sub_species=TimeTree_returned%in%Time_Tree_subspecies$Sub_species_names) %>% 
  left_join(Time_Tree_subspecies,by=c("TimeTree_returned"="Sub_species_names"),suffix=c(".x","_Corresponding")) %>% 
  mutate(Time_Tree_curated=ifelse(is_sub_species,Species_name_Corresponding,TimeTree_returned)) #%>% 
  #select(SampleID,Time_Tree_curated,TimeTree_returned,Species_name_Corresponding,is_sub_species)

# 3.2. Species with only Genus name in TimeTree column: remove it 
# --------------------------------# --------------------------------
genus=c("Artibeus","Equus","Myotis","Nicteris","Rhinolophus")

test = Joined_Metadata %>% 
  subset(TimeTree_returned%in%genus)

# Remove those samples 
Joined_Metadata = Joined_Metadata %>% 
  subset(!TimeTree_returned%in%genus)

# Stats on the new subset 
Ind_per_species <- Joined_Metadata %>%
  group_by(Time_Tree_curated) %>%
  summarise(N_indv=n())
Ind_per_species %>% ggplot(aes(N_indv)) + geom_histogram() + scale_x_log10()
Ind_per_species  %>%
  subset(N_indv>1)  %>%
  pull(N_indv) %>%
  sum() #1557 individuals 

# 3.3. Select species with more than one individual 
# --------------------------------# ---------------
Ind_per_species %>%
  mutate(stat=N_indv>1) %>% pull(stat) %>%  table()

# Prep Subset samples 
Species_list = Ind_per_species %>%
  subset(N_indv>1)  %>%
  pull(Time_Tree_curated)

# Subset samples 
Joined_Metadata <- Joined_Metadata %>% 
  subset(TimeTree_returned%in%Species_list)

Species_from_duplicated_ind = Joined_Metadata  %>% 
  subset(TimeTree_returned%in%unique(Summary_Individual_samples$TimeTree_returned))


# ----------------------
# 4. Add Diet Metadata 
# ----------------------

threshold=70 # threshold to define carnivore, herbivre and insectivores. 
Joined_Metadata=Joined_Metadata %>%
  mutate(Order=Taxonomy_Order) %>%
  mutate(Family=Taxonomy_Family) %>%
  mutate(Diet=ifelse(ET.Diet.PlantO+ET.Diet.Nect+ET.Diet.Fruit+ET.Diet.Seed>threshold,"herbivore","Unknown")) %>%
  mutate(Diet=ifelse(ET.Diet.Vect+ET.Diet.Vend+ET.Diet.Scav+ET.Diet.Vunk+ET.Diet.Vfish>threshold,"carnivore",Diet)) %>%
  mutate(Diet=ifelse(ET.Diet.Inv>threshold,"insectivore",Diet)) %>%
  mutate(Diet=ifelse(Diet=="Unknown" ,"omnivore",Diet))

mismatch =inner_join(Final_metadata,Selected_Samples,by="SampleID") %>%
  subset(is.na(dada2_output_names)) #only keep sampels that have been filtered (other have been excluded because of low number of reads)

# Replace NA by "unknown"
myList <- setNames(lapply(vector("list", ncol(Joined_Metadata)), function(x) x <- "unknown"), names(Joined_Metadata))
Joined_Metadata <- Joined_Metadata %>%  replace_na(myList) 

# Stats 
mytable=function(x){tibble(type=names(table(x)),count=table(x))}

Stats_on_factors=Joined_Metadata %>% 
  select(sample_type,country,Order,Diet,sex,preservative,healthy,LIFE_STAGE,captive_wild,studyID) %>% 
  map(mytable) %>% 
  enframe(name="Factors",value="level") %>% unnest(cols = c(level)) %>% 
  group_by(Factors) %>% 
  mutate(relative_count=round(count/sum(count),2))

a=Joined_Metadata %>% 
  group_by(Order) %>% 
  summarise(n_species_by_order=n_distinct(Time_Tree_curated)) %>% 
  ungroup()

#### SAVING THE OUTPUT METADATA TABLE FOR FURTHER USE 
#####################################################
#write_csv(Joined_Metadata,"Data/Song_Data/Metadata/Subseted_Song_Metadata_withgroup_ID.csv")
write_csv(Joined_Metadata,"Data/Final_Processed_Song_Data/Subseted_Song_Metadata_withgroup_ID.csv")

#Joined_Metadata <- read_csv("Data/Song_Data/Metadata/Subseted_Song_Metadata_withgroup_ID.csv")

length(unique(Joined_Metadata$SampleID))
length(unique(Joined_Metadata$individual_ID))
sum(is.na(Joined_Metadata$individual_ID))
#sum(Ind_per_species$N_indv>1) # 196 species 
#mismatch = Final_metadata %>%
#subset(!SampleID%in%Selected_Samples$SampleID)

#sum(!old_select$SampleID%in%Metadata_online$SampleID)
