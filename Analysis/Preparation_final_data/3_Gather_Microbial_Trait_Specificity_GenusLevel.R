###########################################################################################################################
###########################################################################################################################
#############################################  Specificity and traits at the GENUS LEVEL ##################################
###########################################################################################################################
###########################################################################################################################



rm(list=ls())
set.seed(1)
library(tidyverse)
library(ape)

#########################################
#########################################
####  Specificity at the GENUS LEVEL ####
#########################################
#########################################

# Define thresholds for ASVs
N_sample_min_perASV = 2
N_reads_min_perASV = 100

# Load data 
ASV_specificitySong_raw=as_tibble(read.table("Data/Final_Processed_Song_Data/ASV_specificity_Song_5k.csv",stringsAsFactors = F))
ASV_specificityYoungblut_raw=as_tibble(read.table("Data/Final_Processed_Youngblut_Data/ASV_specificity_5k_YoungBlut.csv",stringsAsFactors = F))

# Filter data 

ASV_specificitySong_raw %>%
  mutate(Abundant=(observed_samples_counts_SONG_morethan2_samples_per_host_5k>N_sample_min_perASV)&(observed_readCounts_SONG_morethan2_samples_per_host_5k>N_reads_min_perASV)) %>% #keep only ASVs observed in at least 3 samples. 
  group_by(Abundant) %>% 
  summarise(nASVs=n(),
            nReads=sum(observed_readCounts_SONG_morethan2_samples_per_host_5k)) %>% 
  mutate(propASV=nASVs/sum(nASVs),
         propReads=nReads/sum(nReads))

ASV_specificitySong_raw %>%
  mutate(Abundant=(observed_samples_counts_SONG_morethan2_samples_per_host_5k>N_sample_min_perASV)&(observed_readCounts_SONG_morethan2_samples_per_host_5k>N_reads_min_perASV)) %>% #keep only ASVs observed in at least 3 samples. 
  mutate(Assignation=!is.na(Genus_SONG_morethan2_samples_per_host_5k)) %>%
  group_by(Assignation,Abundant) %>% 
  summarise(nASVs=n(),
            nReads=sum(observed_readCounts_SONG_morethan2_samples_per_host_5k)) %>% 
  mutate(propASV=nASVs/sum(nASVs)/2,
         propReads=nReads/sum(nReads)/2)

ASV_specificitySong_filt = ASV_specificitySong_raw %>% 
  rename(ASV=seq_SONG_morethan2_samples_per_host_5k)  %>% 
  rename(Genus=Genus_SONG_morethan2_samples_per_host_5k)  %>% 
  rename(Family=Family_SONG_morethan2_samples_per_host_5k)  %>% 
  subset(observed_samples_counts_SONG_morethan2_samples_per_host_5k>N_sample_min_perASV) %>% #keep only ASVs observed in at least 3 samples. 
  subset(observed_readCounts_SONG_morethan2_samples_per_host_5k>N_reads_min_perASV) %>% 
  subset(!is.infinite(sr.obs.z_SONG_morethan2_samples_per_host_5k)) #%>% 
  #left_join(taxo138,by=c("ASV"="ASV_Silva_v138")) %>% 
  #left_join(taxoGTDB,by=c("ASV"="ASV_GTDB"))


# Youngblut 

ASV_specificityYoungblut_raw %>%
  mutate(Abundant = (observed_samples_counts_morethan2_samples_per_host_5k_Youngblut>N_sample_min_perASV)&(observed_readCounts_morethan2_samples_per_host_5k_Youngblut>N_reads_min_perASV)) %>% #keep only ASVs observed in at least 3 samples. 
  group_by(Abundant) %>% 
  summarise(nASVs=n(),
            nReads=sum(observed_readCounts_morethan2_samples_per_host_5k_Youngblut),) %>% 
  mutate(propASV=nASVs/sum(nASVs),
         propReads=nReads/sum(nReads))

ASV_specificityYoungblut_raw %>%
  mutate(Abundant = (observed_samples_counts_morethan2_samples_per_host_5k_Youngblut>N_sample_min_perASV)&(observed_readCounts_morethan2_samples_per_host_5k_Youngblut>N_reads_min_perASV)) %>% #keep only ASVs observed in at least 3 samples. 
  mutate(Assignation=!is.na(Genus_morethan2_samples_per_host_5k_Youngblut)) %>%
  group_by(Assignation,Abundant) %>% 
  summarise(nASVs=n(),
            nReads=sum(observed_readCounts_morethan2_samples_per_host_5k_Youngblut),) %>% 
  mutate(propASV=nASVs/sum(nASVs)/2,
         propReads=nReads/sum(nReads)/2)

ASV_specificityYoungblut_filt = ASV_specificityYoungblut_raw %>% 
  rename(Genus=Genus_morethan2_samples_per_host_5k_Youngblut)  %>% 
  subset(observed_samples_counts_morethan2_samples_per_host_5k_Youngblut>N_sample_min_perASV) %>% #keep only ASVs observed in at least 3 samples. 
  subset(observed_readCounts_morethan2_samples_per_host_5k_Youngblut>N_reads_min_perASV) %>% 
  subset(!is.infinite(sr.obs.z_morethan2_samples_per_host_5k_Youngblut))


# Define variable to summarize
#class_summary=sapply(FUN=class,ASV_specificity)
#numeric_variables=colnames(ASV_specificity)[class_summary%in%c("integer","numeric")]
numeric_variables = c("observed_SR_SONG_morethan2_samples_per_host_5k","observed_samples_counts_SONG_morethan2_samples_per_host_5k",
                      "observed_readCounts_SONG_morethan2_samples_per_host_5k","PD_SONG_morethan2_samples_per_host_5k","sesPD_SONG_morethan2_samples_per_host_5k", 
                      "ASVs_in_Genus_SONG_morethan2_samples_per_host_5k","sr.obs.z_SONG_morethan2_samples_per_host_5k","sr.obs.p_SONG_morethan2_samples_per_host_5k",
                      "ASVs_in_Genus_SONG_morethan2_samples_per_host_5k")


# Compute summary statistics
Genus_Specificity_raw=ASV_specificitySong_filt %>%
 # subset(!is.na(Genus_Silva_v138)) %>%
  group_by(Genus) %>%
  summarise_at(
    .vars = numeric_variables,
    .funs = list(~median(.,na.rm = T),~sum(.,na.rm = T),~mean(.,na.rm = T)))   %>%
  data.frame()
#rownames(Genus_Specificity_raw)=Genus_Specificity_raw$Genus_Silva_v138


# Compute SES PD stats only for ASVs occuring in more than 2 species 
Genus_SpecificityPD_alt = ASV_specificitySong_filt %>%
  subset(!is.na(Genus)) %>%
  subset(observed_SR_SONG_morethan2_samples_per_host_5k>2) %>% 
  select(Genus,sesPD_SONG_morethan2_samples_per_host_5k,sr.obs.z_SONG_morethan2_samples_per_host_5k) %>% 
  group_by(Genus) %>%
  summarise(n=n(),
            sesPD_withoutZero_median = median(sesPD_SONG_morethan2_samples_per_host_5k),
            SR_Z_PDsampling_median=median(sr.obs.z_SONG_morethan2_samples_per_host_5k)) %>%
  data.frame()

Genus_Specificity <- Genus_Specificity_raw %>% 
  left_join(Genus_SpecificityPD_alt)

#Adding YoungBlut estimates 
#class_summary=sapply(FUN=class,ASV_specificityYoungblut)
#numeric_variables=colnames(ASV_specificityYoungblut)[class_summary%in%c("integer","numeric")]
numeric_variables = c("sr.obs_morethan2_samples_per_host_5k_Youngblut","sr.obs.z_morethan2_samples_per_host_5k_Youngblut","sr.obs.p_morethan2_samples_per_host_5k_Youngblut","sesPD_morethan2_samples_per_host_5k_Youngblut" )

# Compute summary statistics
Genus_SpecificityY=ASV_specificityYoungblut_filt %>%
  subset(!is.na(Genus)) %>%
  group_by(Genus) %>%
  summarise_at(
    .vars = numeric_variables,
    .funs = list(~mean(.,na.rm = T),~median(.,na.rm = T),~sum(.,na.rm = T)))   %>%
  data.frame()

rownames(Genus_SpecificityY)=Genus_SpecificityY$Genus

# Full join of the two tables 
Genus_SpecificityBothDB=full_join(Genus_Specificity,Genus_SpecificityY,by="Genus")
Genus_SpecificityBothDB$Genus


#############################
#############################
#   Microbial  traits       #
#############################
#############################

#Load traits
Bergey_Traits_raw <- read_delim("Data/Microbial_trait_Data/Genus_oxygen_spore_Bergey_Mazel.csv",delim = ";")

# define trait categories and reformat database accordingly
non_Anaerobes <- c("Microaerotolerant","Obligately aerobic","Microaerophilic","Aerobicorfacultativelyanaerobic","Aerotolerant", "Facultativeanaerobes","Aerobe")               
StrictAnaerobes <- c("Obligate_anaerobe","obligately anaerobic","Strictlyanaerobic")
Non_StrictAnaerobes <- c("Anaerobe","Obligate_anaerobe","obligately anaerobic","Strictlyanaerobic")
Spore_toNa <- c("spore forming?","Variable","Sporeforming/Variable")

Bergey_Traits  <-  Bergey_Traits_raw %>% 
  mutate(Oxygen_Conservative = ifelse(Oxygen_Manual_Conservative%in%non_Anaerobes,"Non anaerobes",ifelse(Oxygen_Manual_Conservative%in%StrictAnaerobes,"Anaerobes",NA)),
         Oxygen_MajorityRule = ifelse(Oxygen_Manual_MajorityRule%in%non_Anaerobes,"Non anaerobes",ifelse(Oxygen_Manual_MajorityRule%in%Non_StrictAnaerobes,"Anaerobes",NA)),
         Spore_Conservative = ifelse(Spore_Manual_Conservative%in%Spore_toNa,NA,Spore_Manual_Conservative),
         Spore_MajorityRule = ifelse(Spore_Manual_MajorityRule%in%Spore_toNa,NA,Spore_Manual_MajorityRule))   

table(Bergey_Traits$Oxygen_Conservative)
table(Bergey_Traits$Oxygen_MajorityRule)

# BSD 
BSD=read.table("Data/Microbial_trait_Data/Moeller_BSD_ratios.txt",stringsAsFactors = F,header = T)
rownames(BSD)=BSD$Genus

# Add traits to the database
Genus_SpecificityBothDB$BSD=BSD[Genus_SpecificityBothDB$Genus,1]
Genus_SpecificityBothDB_trait <- Genus_SpecificityBothDB %>% 
  left_join(Bergey_Traits)

# Save data 
write.csv(Genus_SpecificityBothDB_trait ,"Data/SubmissionData/Genus_Specificity_Song_Youngblut.csv")






