########################################################
### ------------------------------------------------ ###
###   Download EMP data from mammalian hosts         ###
### -------------------------------------------------###
########################################################


rm(list=ls())
setwd("/Users/fmazel/Desktop/Work/Analyses_en_cours/Species_Specificity/")

###  ----------------------------
###   Load data and funcitons
###  -----------------------------

library(tidyverse)

# Get overall metadata files and sampleIDs
Metadata_online=read_csv("Data/Song_Data/Metadata/inline-supplementary-material-1.csv")
length(unique(Metadata_online$SampleID))

# Get PRJ names 
PRJ_names=readxl::read_xlsx("Data/Song_Data/Metadata/List_studies_ENA_projects.xlsx")


###  ----------------------------
###   Subset samples of mammals 
###  -----------------------------

# Subset data to mammals only 
Metadata_online %>%
  subset(Taxonomy_Class=="Mammalia") %>%
  pull(TimeTree_returned) %>%
  table() %>%
  table()

#select studies 
List_studyID_to_Retrieve=Metadata_online %>%
  subset(Taxonomy_Class=="Mammalia") %>%
  pull(studyID) %>%
  unique() %>%
  as.character 

#List of samples to subset from those studies
List_samples_to_Retrieve=Metadata_online %>%
  subset(Taxonomy_Class=="Mammalia") #%>%
  select(studyID,SampleID)

table(List_samples_to_Retrieve$studyID)
table(Metadata_online$studyID) #2245 Song samples 
length(unique(Metadata_online$SampleID))


dim(List_samples_to_Retrieve)

###  -----------------------------------
###  DELSUC 2014 subsetting and cleaning 
###  -----------------------------------

# Remove contamiated samples from Delsuc et al 2014 
Delsuc_samples=read.csv("Data/Song_Data/Metadata/Delsuc2014_discarded_samples.csv",header = T,stringsAsFactors = F)

Non_Contaminated_Delsuc_samples=Delsuc_samples %>% 
                              subset (!(Unknown>.9999|Soil>.01)) %>%
                              separate(SampleID,sep=".6",into='New_name') %>%
                              mutate(New_name=paste("1056.",New_name,sep="")) %>%
                              mutate(New_name=ifelse(New_name=="1056.Ory.afe.Colch","1056.Ory.afe.Colch.6",New_name)) %>%
                              mutate(New_name=ifelse(New_name=="1056.Myr.tri.Seat.2.03G","1056.Myr.tri.Seat.2.03G06",New_name)) %>%
                              mutate(New_name=ifelse(New_name=="1056.Myr.tri.Seat.1.04J","1056.Myr.tri.Seat.1.04J06",New_name)) %>%
                              pull(New_name)

List_samples_to_Retrieve = List_samples_to_Retrieve %>%
                                subset(!studyID==1056|SampleID%in%Non_Contaminated_Delsuc_samples)

# Adjust sample names 
#s="1056"
#prj=subset(PRJ_names,Qiita_ID==s)$ENA_PRJ
#ftp_adresses=read.table(paste('Data/Song_Data/Accession_files/filereport_read_run_',prj,'_tsv.txt',sep=''),header=T,row.names = NULL,stringsAsFactors = F)
#ftp_adresses = ftp_adresses %>%
#  mutate(sample_title=paste(s,sample_title,sep="."))
#write.table(ftp_adresses,paste('Data/Song_Data/Accession_files/filereport_read_run_',prj,'_tsv.txt',sep=''))


####  -------------------------
#### Download from ENA website 
####  -------------------------


raw_fastq_dir="/Users/fmazel/Data/Mammalian_microbiome_EMP/Raw_fastq/"
done=list.files(raw_fastq_dir)
unknown=c('10376','10925')
new_list=List_studyID_to_Retrieve[!List_studyID_to_Retrieve%in%c(done,unknown)]


new_list="11166"
for (s in new_list) #loop across studies 
  {

  # Create directory for this study
    new_dir=paste(raw_fastq_dir,s,'/',sep="")
    dir.create(new_dir)

  # Load associated ftp adresses 
    prj=subset(PRJ_names,Qiita_ID==s)$ENA_PRJ
    ftp_adresses=read.table(paste('Data/Song_Data/Accession_files/filereport_read_run_',prj,'_tsv.txt',sep=''),header=T,row.names = NULL,stringsAsFactors = F)

  # Subset to samples from EMP list
    ftp_adresses_subset = ftp_adresses %>% 
    subset(sample_title%in%List_samples_to_Retrieve$SampleID)

  #Dowload samples (fastq files from ENA)
    for (sa in 1:dim(ftp_adresses_subset)[1]) #loop across samples within a study
      {
        ftp=as.character(ftp_adresses_subset[sa,"submitted_ftp"])
        sample=strsplit(strsplit(ftp,".fas")[[1]][1],"/")[[1]][6]
        if (paste(sample,".gz",sep="")%in%list.files(new_dir)) {sample=paste(sample,"duplicated",sep="_")} # for multiple sequencing of the same sample 
        download.file(ftp,destfile = paste(new_dir,sample,".gz",sep=""))
      }

  }
   

####  -------------------------
#### Download from MG_RAST 
####  -------------------------

MG_Rast_Sanders=read.csv("Data/Song_Data/Metadata/MG_Rast_IDs_Sanders_Whales.csv",header=T,sep="\t",stringsAsFactors = F)

list_MG_IDs = MG_Rast_Sanders %>%
                subset(type=="Amplicon"&method=="illumina") %>%
                mutate(old_name=name)  %>%
                tidyr::separate(name,sep=".16S.Ilm.V4",into='new_name')  %>%
                mutate(SampleID=paste("10285",new_name,sep='.'))
  

# Subset to samples from EMP list
list_MG_IDs = list_MG_IDs %>% 
  subset(!SampleID%in%List_samples_to_Retrieve$SampleID)

dir_Sanders='/Users/fmazel/Data/Mammalian_microbiome_EMP/Raw_fastq/10285/'

for (i in 1:dim(list_MG_IDs)[1])
{
  MG_adress=paste("http://api.metagenomics.anl.gov/1/download/",list_MG_IDs$MG.RAST_ID[i],"?file=050.1",sep="")
  destdir=paste(dir_Sanders,list_MG_IDs$new_name[i],".fastq",sep="")
  download.file(MG_adress,destfile = destdir)
}

###  ------------------------------------------------
###   Remove reads from very deeply sequenced samples -- organize the outputs per studies 
###  -------------------------------------------------

# BASH script  

ls /Users/fmazel/Data/Mammalian_microbiome_EMP/Raw_fastq/ | cat | while read studies 
do
mkdir /Users/fmazel/Data/Mammalian_microbiome_EMP/Non_Combined_Rarefied_fastq/"$studies"
done

fmazel/seqtk/seqtk 
ls /Users/fmazel/Data/Mammalian_microbiome_EMP/Raw_fastq/ | cat | while read studies 
do
ls /Users/fmazel/Data/Mammalian_microbiome_EMP/Raw_fastq/"$studies" | cat | while read line 
do
echo "$line"
fmazel/seqtk/seqtk sample -s100 -s100 /Users/fmazel/Data/Mammalian_microbiome_EMP/Raw_fastq/"$studies"/"$line" 20000 > /Users/fmazel/Data/Mammalian_microbiome_EMP/Non_Combined_Rarefied_fastq/"$studies"/"$line" ;
done

done 

###  --------------------------------------------
###   Reverse Complement study from Sanders et al 
###  -------------------------------------------

ls /Users/fmazel/Data/Mammalian_microbiome_EMP/Non_Combined_Rarefied_fastq/10285 | cat | while read line 
do
echo "$line"
/Users/fmazel/seqtk/seqtk seq -r /Users/fmazel/Data/Mammalian_microbiome_EMP/Non_Combined_Rarefied_fastq/10285/"$line"  > /Users/fmazel/Data/Mammalian_microbiome_EMP/Non_Combined_Rarefied_fastq/10285_RC/"$line"
done





