# Gather metadata 


#List of ENA IDs
#OTUtable=read.table("/Users/fmazel/Documents/CoPhylogeny/Dada2_ASVs/ASVs_counts.txt")
#samplesENA=colnames(OTUtable)

#List of ENA IDs
ENA_output=read.table("Data/Youngblut_Data/PRJEB29403.txt",header=T,stringsAsFactors = F)
ENA_output$sampleIDs=sapply(ENA_output$sra_ftp,FUN=function(x){strsplit(strsplit(strsplit(x,";")[[1]][1],"/")[[1]][6],".R1")[[1]][1]})
ENA_output$ERR_name=sapply(ENA_output$sra_ftp,FUN=function(x){strsplit(strsplit(x,";")[[1]][1],"/")[[1]][5]})
table(ENA_output$sampleIDs)[table(ENA_output$sampleIDs)>1]

#remove non singletons
toremove=c("298_Horsfields_Bronze_Cuckoo","398_Eurasian_Reed_Warbler","F300_Asp")
ENA_output=ENA_output[!ENA_output$sampleIDs%in%toremove,]
ENA_output=ENA_output[!ENA_output$ERR_name=="ERR2860829",] #remove one sample from Kangourou 

rownames(ENA_output)=ENA_output$sampleIDs

#Metadata
library(openxlsx)
Metadata=read.xlsx("Data/Youngblut_Data/41467_2019_10191_MOESM5_ESM.xlsx")
Metadata$sample_ID_2=substring(gsub(x=Metadata$sample_ID,pattern="\\.",replacement = '_',),2)


#remove the 'F' in front of some samples 
ENA_output$sampleIDs_V2=ENA_output$sampleIDs
ENA_output$sampleIDs_V2[substring(ENA_output$sampleIDs,1,1)=="F"]=substring(ENA_output$sampleIDs[substring(ENA_output$sampleIDs,1,1)=="F"],2)


rownames(Metadata)=Metadata$sample_ID_2
rownames(ENA_output)=ENA_output$sampleIDs_V2
Metadata$ENA_Id=ENA_output[rownames(Metadata),"ERR_name"]

write.table(Metadata,file="Data/Youngblut_Data/Updated_Metadata.txt")


