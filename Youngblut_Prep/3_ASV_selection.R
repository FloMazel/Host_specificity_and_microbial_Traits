####################################################
######   READS processing Youngblut et al.    ######
####################################################


rm(list=ls())
library(dada2)
library(tidyverse)

path <- "/Users/fmazel/Documents/CoPhylogeny/RawData/" 
pathFilt <- "/Users/fmazel/Documents/CoPhylogeny/Dada2_Filtered/" 

fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

####plot quality profiles for forward and reverse reads ####
plotQualityProfile(fnFs[1:20])
plotQualityProfile(fnRs[1:20])

####trim & filter####
filtFs <- file.path(pathFilt, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(pathFilt, "filtered", paste0(sample.names, "_R_filt.fastq"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,matchIDs=F,
                     compress=TRUE, multithread=TRUE) #the more aggressive you are with truncation, the more reads will be retained, but you still need to merge, so consider how much read overlap remains after truncation

percent_retained <- as.data.frame(out)
percent_retained$percentage <- (percent_retained[,2]/percent_retained[,1])*100
hist(percent_retained$percentage)

####dereplication####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

####learn error rates####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

####sample inference####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

####merge paired reads####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE) #a "subscript out of bounds" error here may indicate that you aren't merging any reads in one or more samples. you can remove samples with low counts from the workflow before the filterAndTrim step (a few steps back)

####construct sequence table####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#look at sequence length distribution (253 on average)
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
plot(x=length.histogram[,1], y=length.histogram[,2]) #view length distribution plot

###chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=T, verbose=T) #
dim(seqtab.nochim)

# though we only lost 21 sequences, we don't know if they held a lot in terms of abundance, this is one quick way to look at that
sum(seqtab.nochim)/sum(seqtab) # 

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {asv_headers[i] <- paste(">ASV", i, sep="_")}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste(path,"/ASVs.fa",sep=""))

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste(path,"/ASVs_counts.txt",sep=""), sep="\t", quote=F)

#assigning taxonomy
taxa <- assignTaxonomy(seqs=colnames(seqtab.nochim), refFasta="/Users/fmazel/Data/DataBase/silva_nr_v132_train_set.fa", multithread=T, tryRC=T)
r=data.frame(taxa)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, paste(path,"/ASVs_taxonomy.txt",sep=""), sep="\t", quote=F)

# ----------------------
#  ASVs filtering 
# ----------------------
# We only keep bacterial ASVs with a known phylum and exclude mitochondria, chloroplasts and rare ASVs (less than 10 reads or only one sample)
asv_table <- read.table("Data/Youngblut_Data/ASVs_counts.txt")
ASV_taxonomy <- read.table("Data/Youngblut_Data/ASVs_taxonomy.txt")
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

# ----------------------
#  Export Filtered data
# ----------------------

path_to_final_results = "Data/Final_Processed_Youngblut_Data/"

# count table
asv_tab_filt = asv_table[ASV_taxonomy_filt$ASV,]
write.table(asv_table, paste(path_to_final_results,"/ASVs_counts_filt.txt",sep=""), sep="\t", quote=F)

# tax table:
write.table(ASV_taxonomy_filt, paste(path_to_final_results,"/ASVs_taxonomy_filt.txt",sep=""), sep="\t", quote=F)

# Rarefied table
depth=apply(asv_tab_filt,2,sum) #check depth distribution 
asv_tab_rar <- rrarefy(t(asv_tab_filt),4900)
depth_rar=apply(asv_tab_rar,1,sum) #check depth distribution 
asv_tab_rar <- asv_tab_rar[depth_rar==4900,] #removing sample with < than 5000 reads

ASV_prev=apply(asv_tab_rar,2,sum)
asv_tab_rar <- asv_tab_rar[,ASV_prev>0] #removing ASV with no occurence (becauise of rarefaction, basically)
dim(asv_tab_rar)

write.table(asv_tab_rar, paste(path_to_final_results,"/ASVs_counts_filt_rar4900.txt",sep=""), sep="\t", quote=F)

