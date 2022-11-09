#define plot features 

# 1 pt = 0.35mm
# Recommended max size=7 , min size =5 in pt 
P1=5;P2=6;P3=7
MyTheme=theme_classic() + 
  theme(plot.title = element_text(size=P3),
        axis.text=element_text(size=P2),
        axis.text.x =element_text(size=P2,angle=45,hjust=1),
        axis.title=element_text(size=P2),
        
        legend.text=element_text(size=P2),
        legend.title=element_text(size=P2),
        
        panel.border = element_blank(),
        
        axis.ticks = element_line(size = 1*.35),
        axis.ticks.length = unit(.5, "mm"),
        
        plot.caption = element_text(hjust = 0), 
        plot.title.position = "plot")

mytable=function(x){tibble(type=names(table(x)),count=table(x))}


re_shape_name=function(x){
  splitted_name=strsplit(x,".gz")
  initial_name=strsplit(x,"_F_filt.fastq")[[1]][1]
  tibble(initial_fastq_name =initial_name,dada2_output_names=x,SampleID=splitted_name[[1]][1],Other=splitted_name[[1]][2])
}



# Function that change the name of the entry of the supp tsable presenting the results of the permanova 
changeSampleName = Vectorize(function(x){
  if (x=="Samples_from_species_with_5_samples") {y="5 individuals per species"}
  if (x=="Samples_from_species_with_more_than_1_sample") {y="All individuals"}
  return(y)
})
changeBetaName = Vectorize(function(x){
  if (x=="Bray") {y="Bray-curtis"}
  else if (x=="Jac") {y="Jaccard"}
  else if (x=="Jac_TT") {y="True_turnover_Jaccard"}
  return(y)
})

# This function load fastq files in R and extract reads info (length, counts)
import_analyze_fastq=function(study_name,fastq_name,raw_fastq_dir){
  
  fastq_dir = paste(raw_fastq_dir,study_name,sep="")
  fastq_file=paste(fastq_dir,fastq_name,sep="/")
  
  imported_fastq=readFastq(fastq_file)
  fastq=imported_fastq@sread@ranges@width
  
  results=tibble(
    studyID=study_name,
    sample_name=fastq_name,
    minLength=min(fastq),
    maxLength=max(fastq),
    medianLength=median(fastq),
    meanLength=mean(fastq),
    read_number=length(imported_fastq@sread)
  )
  
  return(results)
}


# MAke a list of all sampels and their associated study based on raw fastq file dowloaded 
makeSampleList=function(raw_fastq_dir){
  res=list()
  list_studies=list.files(raw_fastq_dir)
  for (s in list_studies) {
    fastq_dir=paste(raw_fastq_dir,s,sep="")
    res[[s]]=tibble(
      studyID=s,
      sampleID=list.files(fastq_dir),
    )
    
    res[[s]]= res[[s]] %>% 
      mutate(file_path=paste(raw_fastq_dir,studyID,'/',sampleID,sep=''))
    
  }
  
  All_samples=do.call(rbind,res)
  return(All_samples)
}


# Define marginal histogram
marginal_distribution <- function(x, var, group) {
  ggplot(x, aes_string(x = var, fill = group)) +
    geom_histogram(bins = 20, alpha = 0.4, position = "identity") +
    # geom_density(alpha = 0.4, size = 0.1) +
    guides(fill = FALSE) +
    theme_void() +
    scale_x_continuous(trans='log')+
    theme(plot.margin = margin())
}

# Function to wrap adonis and run it on multiple datasets
wrapper_adonis=function(i,Disp_results,Betas,Sample_list)
{
  print(i)
  
  # Retrieve parameters 
  # --------------------
  Phyl_resolution=Disp_results$Phylo_resolution[i]
  Beta_metric=Disp_results$Metrics[i]
  Sample_set=Disp_results$Sample_set[i]
  Data_set=Disp_results$Data_set[i]
  Metadata=Sample_list[[Data_set]][[Sample_set]]
  Strata=Disp_results$Strata[i]
  
  # Prepare data 
  # --------------------
  Beta_subset=Betas[[Data_set]][[Phyl_resolution]][Beta_metric,as.character(Metadata$SampleID),as.character(Metadata$SampleID)]
  #Host_species=Metadata$Time_Tree_curated
  #print(head(Host_species))
  
  
  # Run the model and export output
  # --------------------
  if (Strata=="None"){output=adonis2(Beta_subset~Time_Tree_curated,data=Metadata) } else {
    Strata_list=pull(Metadata[,Strata])
    output=adonis2(Beta_subset~Time_Tree_curated,data=Metadata,strata=Strata_list) }
  
  #res= c(output$aov.tab[1,"F.Model"],output$aov.tab[1,"Pr(>F)"],output$aov.tab[1,"R2"],output$aov.tab[1,"Df"]+1) #for adonis
  res= c(output$F[1],output[["Pr(>F)"]][1],output$R2[1],output$Df[1]+1)#for adonis2
  
  return(res)
}

# function copied from https://kieranhealy.org/blog/archives/2018/11/06/spreading-multiple-values/

myspread <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

# Function to wrap betadispers across datasets 
wrapper_betadisper=function(i,Disp_results,Betas,Sample_list)
{
  print(i)
  
  # Retrieve parameters 
  # --------------------
  Phyl_resolution=Disp_results$Phylo_resolution[i]
  Beta_metric=Disp_results$Metrics[i]
  Sample_set=Disp_results$Sample_set[i]
  Data_set=Disp_results$Data_set[i]
  Metadata=Sample_list[[Data_set]][[Sample_set]]
  Strata=Disp_results$Strata[i]
  
  # Prepare data 
  # --------------------
  Beta_subset=(Betas[[Data_set]][[Phyl_resolution]][Beta_metric,as.character(Metadata$SampleID),as.character(Metadata$SampleID)])
  Host_species=Metadata$Time_Tree_curated
  
  
  # Run the model and export output
  # --------------------
  output=betadisper(as.dist(Beta_subset),Host_species)
  if (Strata=="None"){perms=permutest(output) } else {
    Strata_list=Metadata[,Strata]
    perms=permutest(output,srata=Strata_list)}
  
  res= c(perms$tab[1,"F"],perms$tab[1,"Pr(>F)"],perms$tab[1,"Df"]+1)
  return(res)
}


#funciton to comute specifity 
# arguments: 

specificity=function(subset_name,prevASV_table,Mammals_samples,sample_to_host_species,HostTree,Bacterial_Taxo,runs=100,null.model="richness")
{
  
  # prep data 
  prevASV_table=prevASV_table[,as.character(Mammals_samples$SampleID)] #keep only coorect samples 
  prevASV_table=prevASV_table[apply(prevASV_table,1,sum)>0,] #remove ASVs than are not present in the reduced subset of samples 
  sample_to_host_species=sample_to_host_species[colnames(prevASV_table),unique(as.character(Mammals_samples$Time_Tree_curated))]   # common ASV list
  Sample_distribution_host_species=apply(sample_to_host_species,2,sum) 
  
  # Compute each ASV relative prevalence per host species
  prevASV_table_PA=prevASV_table
  prevASV_table_PA[prevASV_table_PA>0]=1
  Prevalence_Table_Host_species=prevASV_table_PA%*%as.matrix(sample_to_host_species) #compute number of apperance of each ASV in each host species
  Rel_Prevalence_Table_Host_species=t(t(Prevalence_Table_Host_species)/Sample_distribution_host_species[colnames(Prevalence_Table_Host_species)]) # Normalisation step: divide number of apperance of each ASV in each host species by the number of samples per species 
  
  # Compute each ASV relative abundance per species
  Abundance_Table_Host_species=prevASV_table%*%as.matrix(sample_to_host_species) #compute total # reads per ASV per host species
  Tot_reads_per_ASV=apply(Abundance_Table_Host_species,1,sum)
  rel_Abundance_Table_Host_species=Abundance_Table_Host_species/apply(Abundance_Table_Host_species,1,sum) # Normalisation: divide by total number of reads per ASVs
  
  # Compute ASV specificty metrics
  ASV_specificity=tibble(seq=row.names(Prevalence_Table_Host_species),
                         observed_SR=apply(Rel_Prevalence_Table_Host_species>0,1,sum),
                         observed_samples_counts=apply(Prevalence_Table_Host_species,1,sum),
                         observed_readCounts=apply(Abundance_Table_Host_species,1,sum),
                         PD=pd.query(HostTree, Rel_Prevalence_Table_Host_species>0, standardize = FALSE),
                         sesPD=pd.query(HostTree, Rel_Prevalence_Table_Host_species>0, standardize = TRUE),
                         simpson_prev=diversity(Rel_Prevalence_Table_Host_species,"simpson"),
                         shannon_prev=diversity(Rel_Prevalence_Table_Host_species,"shannon"),
                         simpson_ab=diversity(rel_Abundance_Table_Host_species,"simpson"),
                         shannon_ab=diversity(rel_Abundance_Table_Host_species,"shannon")
  )
  
  ASV_specificity$logSR=log(ASV_specificity$observed_SR)
  ASV_specificity$pielou_prev=ASV_specificity$shannon_prev/ASV_specificity$logSR
  ASV_specificity$pielou_ab=ASV_specificity$shannon_ab/ASV_specificity$logSR
  sesSR=ses_SR(prevASV_table,sample_to_host_species,runs=runs,null.model=null.model) #Null model for SR
  
  ASV_specificity=bind_cols(ASV_specificity,sesSR)
  
  # Add taxonomy information 
  
  ASV_specificity$Genus=Bacterial_Taxo[ASV_specificity$seq,"Genus"]
  ASV_specificity$Family=Bacterial_Taxo[ASV_specificity$seq,"Family"]
  ASV_specificity$Order=Bacterial_Taxo[ASV_specificity$seq,"Order"]
  ASV_specificity$Class=Bacterial_Taxo[ASV_specificity$seq,"Class"]
  ASV_specificity$Phylum=Bacterial_Taxo[ASV_specificity$seq,"Phylum"]
  
  # Number of ASV in a particular taxo rank 
  ASV_specificity$ASVs_in_Genus=table(ASV_specificity$Genus)[ASV_specificity$Genus]
  ASV_specificity$ASVs_in_Family=table(ASV_specificity$Family)[ASV_specificity$Family]
  ASV_specificity$ASVs_in_Order=table(ASV_specificity$Order)[ASV_specificity$Order]
  
  rownames(ASV_specificity)=ASV_specificity$seq
  colnames(ASV_specificity)=paste(colnames(ASV_specificity),subset_name,sep="_")
  
  return(ASV_specificity)
}


#function to compute SR
SR=function(OTU_table,sample_by_host_taxonomy_table)
{
  Prevalence_Table_Host_species=OTU_table%*%as.matrix(sample_by_host_taxonomy_table) #compute number of apperance of each ASV in each host species
  observed_SR=apply(Prevalence_Table_Host_species>0,1,sum)
  return(observed_SR)
}

#OTU_table=OTU_table_nonrare_ASV_MMF
#sample_by_host_taxonomy_table=sample_to_host_species
#runs=10

# Function to compute specificity and null models for richness 
ses_SR=function(OTU_table,sample_by_host_taxonomy_table,runs=100,null.model)
{
  
  #convert to presence/absence 
  OTU_tablePA=OTU_table
  OTU_tablePA[OTU_tablePA>0]=1
  
  sr.obs=SR(OTU_tablePA,sample_by_host_taxonomy_table)
  #
  randomizeMatrix(OTU_tablePA,null.model = null.model)
  sr.rand = t(replicate(runs, as.vector(SR(randomizeMatrix(OTU_tablePA,null.model = null.model), sample_by_host_taxonomy_table))))
  
  sr.rand.mean <- apply(X = sr.rand, MARGIN = 2, FUN = mean, 
                        na.rm = TRUE)
  sr.rand.sd <- apply(X = sr.rand, MARGIN = 2, FUN = sd, 
                      na.rm = TRUE)
  sr.obs.z <- (sr.obs - sr.rand.mean)/sr.rand.sd
  sr.obs.rank <- apply(X = rbind(sr.obs, sr.rand), MARGIN = 2, 
                       FUN = rank)[1, ]
  sr.obs.rank <- ifelse(is.na(sr.rand.mean), NA, sr.obs.rank)
  results=data.frame(sr.obs, sr.rand.mean,  sr.rand.sd, sr.obs.rank, sr.obs.z, sr.obs.p = sr.obs.rank/(runs + 1), runs = runs, row.names = row.names(OTU_table))
  return(results)
  
  
}

#a=ses_SR(OTU_table_nonrare_ASV_MMF,sample_by_host_taxonomy_table,runs=100)
