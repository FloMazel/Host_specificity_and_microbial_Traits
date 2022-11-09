########################################
# DADA 2 pipeline for the EMP studies  #
########################################

# ----------------------------
# Load packages and functions 
# ----------------------------

rm(list=ls())
library(dada2)
library(tidyverse)
library(vegan)
source("Scripts/Final/0_Utility_functions.R")

# ----------------------------
# Define raw data directory 
# ----------------------------

n=4 #Number of available core to use (parrallel processing of multiple dada2 steps)
path_to_exported_dada2_output="/Users/fmazel/Data/Mammalian_microbiome_EMP/Dada2_output_separated_studies/"
list_studies=list.files("/Users/fmazel/Data/Mammalian_microbiome_EMP/Raw_fastq/")
pathraw <- "/Users/fmazel/Data/Mammalian_microbiome_EMP/Non_Combined_Rarefied_fastq/"
pathFilt <- "/Users/fmazel/Data/Mammalian_microbiome_EMP/Non_Combined_Dada2_Filtered/"

# ----------------------------
#        Infer ASVs 
# ----------------------------

# We loop over studeis bc the seqeuncign error might differ between studies. 

for (studies in rev(list_studies))
{
  #path_to_exported_dada2_output_S=paste(path_to_exported_dada2_output,studies,sep="")
  #dir.create(path_to_exported_dada2_output_S)
  #path_to_exported_dada2_output_S=paste(path_to_exported_dada2_output_S,"/",sep="")
  pathFilt_S= paste("/Users/fmazel/Data/Mammalian_microbiome_EMP/Non_Combined_Dada2_Filtered/",studies,sep="")
  #rarefied sequences
  raw_fastq_dir=paste(pathraw,studies,'/',sep='')
  fnFs=list.files(raw_fastq_dir,full.names = T)
  sample.names <- list.files(raw_fastq_dir)
  
  ####plot quality profiles 
  plotQualityProfile(fnFs[1:5])
  
  ####trim & filter####
  filtFs <- file.path(pathFilt_S, "filtered", paste0(sample.names, "_F_filt.fastq"))
  
  out <- filterAndTrim(fnFs, filtFs, truncLen=c(90),
                       maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,matchIDs=F,
                       compress=TRUE, multithread=n) #the more aggressive you are with truncation, the more reads will be retained, but you still need to merge, so consider how much read overlap remains after truncation
  
  percent_retained <- as.data.frame(out)
  percent_retained$percentage <- (percent_retained[,2]/percent_retained[,1])*100
  hist(percent_retained$percentage)
  percent_retained $sample.names=rownames(percent_retained)
  
  #remove samples with less than 5000 reads
  filtered_sample.names = percent_retained %>% 
    subset (reads.out>4999) %>%
    pull(sample.names)
  
  filtFs <- list.files(paste(pathFilt_S,'filtered',sep='/'),full.names = T)
  filtered_sample.names <- list.files(paste(pathFilt_S,'filtered',sep='/'))
  
  ####dereplication####
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  
  # Name the derep-class objects by the sample names
  names(derepFs) <- filtered_sample.names
  
  # Check for problematic samples 
  minQ <- sapply(derepFs, function(x) min(as.vector(x$quals)))
  
  if (sum(minQ < 0)>0)
  {
    to_remove=which(minQ < 0)
    #Only keep the good ones 
    filtered_derepFs=derepFs[names(derepFs)[-to_remove]]
    filtered_filtFs=filtFs[-to_remove]
  } else {
    filtered_derepFs=derepFs
    filtered_filtFs=filtFs
  }
  
  
  ####learn error rates####
  errF <- learnErrors(filtered_filtFs, multithread=n)
  plotErrors(errF, nominalQ=TRUE)
  
  ####sample inference####
  dadaFs <- dada(filtered_derepFs, err=errF, multithread=n)
  
  saveRDS(dadaFs,file=paste(path_to_exported_dada2_output,"DadaRawOutput/Dada_output",studies,".RDS",sep=""))
  
}

# ------------------------
# Combine raw dada outputs 
# ------------------------

data_raw_dir=paste(path_to_exported_dada2_output,"DadaRawOutput/",sep="")
dadaFs_combined=unlist(sapply(list.files(data_raw_dir,full.names = T),readRDS),recursive = FALSE)

sp_names <- tibble(raw_sample_names = names(dadaFs_combined)) %>% 
  separate(raw_sample_names,sep=".RDS.", into = c("A","Fastq_file_name"),remove = F)

# ----------------------------
#  Construct sequence table
# ----------------------------

seqtab <- makeSequenceTable(dadaFs_combined)
dim(seqtab_combined)

#look at sequence length distribution (253 on average)
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
plot(x=length.histogram[,1], y=length.histogram[,2]) #view length distribution plot

# ----------------
#  Chimera removal
# ----------------

seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=T, verbose=T) #
rownames(seqtab.nochim) <- sp_names$Fastq_file_name

dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # 

# --------------------
# Change ASV names
# --------------------

#giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

# keep unique names instead
for (i in 1:dim(seqtab.nochim)[2]) {asv_headers[i] <- paste(">ASV", i, sep="_")}

SP_samples=read_csv("Data/Song_Data/Metadata/Subseted_Song_Metadata_withgroup_ID.csv")

# ----------------------------
# Save raw outputs 
# ----------------------------

path_to_final_results = "Data/Final_Processed_Song_Data/"

#making and writing out a fasta of our final ASV seqs

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste(path_to_final_results,"/ASVs.fa",sep=""))

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste(path_to_final_results,"/ASVs_counts.txt",sep=""), sep="\t", quote=F)

# ----------------------
#  Assigning taxonomy
# ----------------------
taxa <- assignTaxonomy(seqs=colnames(seqtab.nochim), refFasta="/Users/fmazel/Data/SILVA/dada2_formated/silva_nr_v132_train_set.fa", multithread=T, tryRC=T)
r=data.frame(taxa)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, paste(path_to_final_results,"/ASVs_taxonomy.txt",sep=""), sep="\t", quote=F)

taxo132 <- read.table("Data/Final_Processed_Song_Data/ASVs_taxonomy.txt",row.names = NULL)

taxa138 <- assignTaxonomy(seqs=fastaASV, refFasta="/Users/fmazel/Data/SILVA/dada2_formated/silva_nr99_v138.1_train_set.fa.gz", multithread=4, tryRC=T)
row.names(taxa138) <- taxo132$row.names
write.table(taxa138, "Data/Final_Processed_Song_Data/ASVs_taxonomy_SILVA138.txt", sep="\t", quote=F)

taxaGTDB <- assignTaxonomy(seqs=fastaASV, refFasta="/Users/fmazel/Data/dada2_formated_databases/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz", multithread=4, tryRC=T)
row.names(taxaGTDB) <- taxo132$row.names
write.table(taxaGTDB, "Data/Final_Processed_Song_Data/ASVs_taxonomy_GTDB.txt", sep="\t", quote=F)


# ----------------------
#  ASVs filtering 
# ----------------------
# We only keep bacterial ASVs with a known phylum and exclude mitochondria, chloroplasts  and rare ASVs (less than 10 reads or only one sample)

asv_table <- read.table("Data/Final_Processed_Song_Data/ASVs_counts.txt",check.names = F)
ASV_taxonomy <- read.table("Data/Final_Processed_Song_Data/ASVs_taxonomy.txt")
ASV_taxonomy$ASV=rownames(ASV_taxonomy)

counts_stats <- tibble(ASV=rownames(asv_table),
                       prevalence = apply(asv_table>0,1,sum),
                       abundance = apply(asv_table,1,sum))

ASV_taxonomy <- ASV_taxonomy %>% 
  left_join(counts_stats)
# Some stats on reads 
ASV_taxonomy %>% 
  group_by(Kingdom) %>% 
  summarise(Tot_ab = sum(abundance)) %>% 
  mutate(rel_ab = Tot_ab/sum(Tot_ab))

ASV_taxonomy %>% 
  group_by(Order) %>% 
  summarise(Tot_ab = sum(abundance)) %>%
  mutate(rel_ab = Tot_ab/sum(Tot_ab))  %>% 
  subset(Order=="Chloroplast")

ASV_taxonomy %>% 
  group_by(Family) %>% 
  summarise(Tot_ab = sum(abundance)) %>% 
  mutate(rel_ab = Tot_ab/sum(Tot_ab))  %>% 
  subset(Family=="Mitochondria")

sum(ASV_taxonomy$prevalence==1)
sum(ASV_taxonomy$abundance<10)

ASV_taxonomy %>% 
  group_by(Family) %>% 
  summarise(Tot_ab = sum(abundance)) %>% 
  mutate(rel_ab = Tot_ab/sum(Tot_ab))  %>% 
  subset(Family=="Mitochondria")

ASV_taxonomy_filt <- ASV_taxonomy %>% 
  subset(Kingdom == "Bacteria") %>% 
  subset(!Order=="Chloroplast") %>% 
  subset(!Family=="Mitochondria") # %>% 
#subset(prevalence>1)  %>% 
#subset(abundance>10)  

ASV_taxonomy_filt_NoRare <- ASV_taxonomy %>% 
  subset(Kingdom == "Bacteria") %>% 
  subset(!Order=="Chloroplast") %>% 
  subset(!Family=="Mitochondria")  %>% 
  subset(prevalence>1)  
#subset(abundance>10)  

# ----------------------
#  Export Filtered data
# ----------------------
path_to_final_results = "Data/Final_Processed_Song_Data/"

# count table
asv_tab_filt = asv_table[ASV_taxonomy_filt$ASV,]
write.table(asv_tab_filt, paste(path_to_final_results,"/ASVs_counts_filt.txt",sep=""), sep="\t", quote=F)

asv_tab_filtNorare = asv_table[ASV_taxonomy_filt_NoRare$ASV,]
write.table(asv_tab_filtNorare, paste(path_to_final_results,"/ASVs_counts_filt_noRare.txt",sep=""), sep="\t", quote=F)

# tax table:
write.table(ASV_taxonomy_filt, paste(path_to_final_results,"/ASVs_taxonomy_filt.txt",sep=""), sep="\t", quote=F)

# Rarefied table
depth=apply(asv_tab_filtNorare,2,sum) #check depth distribution 
asv_tab_rar <- rrarefy(t(asv_tab_filtNorare),5000)
depth_rar=apply(asv_tab_rar,1,sum) #check depth distribution 
asv_tab_rar <- asv_tab_rar[depth_rar==5000,] #removing sample with < than 5000 reads

ASV_prev=apply(asv_tab_rar,2,sum)
asv_tab_rar <- asv_tab_rar[,ASV_prev>0] #removing ASV with no occurence (becauise of rarefaction, basically)
dim(asv_tab_rar)

write.table(t(asv_tab_rar), paste(path_to_final_results,"/ASVs_counts_filt_rar5000_NoRare.txt",sep=""), sep="\t", quote=F)

# Rarefied table Without Rare 
depth=apply(asv_tab_filt,2,sum) #check depth distribution 
asv_tab_rar <- rrarefy(t(asv_tab_filt),5000)
depth_rar=apply(asv_tab_rar,1,sum) #check depth distribution 
asv_tab_rar <- asv_tab_rar[depth_rar==5000,] #removing sample with < than 5000 reads

ASV_prev=apply(asv_tab_rar,2,sum)
asv_tab_rar <- asv_tab_rar[,ASV_prev>0] #removing ASV with no occurence (becauise of rarefaction, basically)
dim(asv_tab_rar)

write.table(t(asv_tab_rar), paste(path_to_final_results,"/ASVs_counts_filt_rar5000.txt",sep=""), sep="\t", quote=F)



