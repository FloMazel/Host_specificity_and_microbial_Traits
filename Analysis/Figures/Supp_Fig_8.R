## ---------------------------------------- ##
# LINK TO CONSERVATIVE OXYGEN ASSIGNATION
## -----------------------------------------

outputFig <- "Submission/Mol_Ecol/V2/Supp_data/Supp_Figure_8.pdf"


TaxoSpecificitynameS = "SEA data: Taxonomic Specificity (Z-score)"
TaxoSpecificitynameY = "YEA data: Taxonomic Specificity (Z-score)"


# Load data with 

#############################
######    Results      ######
#############################

# --------------------#
# -----------------
# Select the datasets (Y Vs Song)  the metric (host range vs host phylo range)
# -----------------
# --------------------#
dataGenus = Genus_Specificity %>% 
  rename(zscore=sr.obs.z_SONG_morethan2_samples_per_host_5k_median,
         zscore_altforPD=SR_Z_PDsampling_median,
         zscorePD=sesPD_withoutZero_median) %>% 
  dplyr::select(Genus,Genus_Phylo,n, Phylum_simplified,Phylum, zscore,zscorePD,zscore_altforPD,BSD,
         Oxygen_Conservative,Oxygen_MajorityRule,Spore_Conservative,Spore_MajorityRule,ASVs_in_Genus_SONG_morethan2_samples_per_host_5k_sum)

# Sample size
Data_ox <- dataGenus %>% 
  mutate(Oxygen_trait = Oxygen_Conservative) %>%  
  dplyr::select(Oxygen_trait,zscore,Genus,Genus_Phylo,Phylum,Phylum_simplified) %>% 
  drop_na()
dim(Data_ox)[1]



# PGLS Model
nTree=length(myGenusTrees)
Pgls_outupt = tibble(F_val=rep(NA,nTree),p=rep(NA,nTree),lambda=rep(NA,nTree))

for (i in 1:nTree) {
  dataCaper<- comparative.data(myGenusTrees[[i]],  Data_ox, Genus, vcv=TRUE, vcv.dim=3)
  PGLS <- pgls(zscore ~ Oxygen_trait, data = dataCaper, lambda='ML')
  aov_PGLS =  anova.pgls(PGLS) 
  
  Pgls_outupt[i,"F_val"] = aov_PGLS$`F value`[1]
  Pgls_outupt[i,"p"] = aov_PGLS$`Pr(>F)`[1]
  Pgls_outupt[i,"lambda"] = PGLS$param["lambda"]
}

median_PGLS_output = apply(Pgls_outupt,2,median)
median_PGLS_output 
sd_PGLS_output = apply(Pgls_outupt,2,sd)
sd_PGLS_output


if (median_PGLS_output["p"]<.001) {kp="<.001"} else {kp=round(median_PGLS_output["p"],3)}


PGLS_model_text = paste("F value=",round(median_PGLS_output["F_val"],2),
                        "\np=",kp,
                        "; n=",dim( Data_ox)[1],
                        sep="")

OxygenPlotS = Data_ox %>% 
  subset(!is.na(Oxygen_trait)) %>% 
  ggplot(aes(y=-zscore,
             x=Oxygen_trait,na.rm = T))+ 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(col=Phylum_simplified),alpha=alphaJitter,size=sizeJitter) + scale_color_manual(values=cols)  +
  ylab(TaxoSpecificitynameS) + xlab("Oxygen requirement")+
  annotate(geom="text", x=1.5,y=Ystats , label=PGLS_model_text ,size=P2*0.35) +
  ylim(c(YlimDown,YlimUp )) +
  MyTheme + theme(legend.position="none")+ ggtitle("SEA dataset")

OxygenPlotS




## 3. Supp plot: Younggblut AND HOST range
## ----------------------------------## 
dataGenus = Genus_Specificity %>% 
  rename(zscore=sr.obs.z_morethan2_samples_per_host_5k_Youngblut_median,
         zscorePD=sesPD_morethan2_samples_per_host_5k_Youngblut_median) %>% 
  mutate(zscore_altforPD=zscore) %>% 
  dplyr::select(Genus,Genus_Phylo,n, Phylum_simplified,Phylum, zscore,zscorePD,zscore_altforPD,BSD,
         Oxygen_Conservative,Oxygen_MajorityRule,Spore_Conservative,Spore_MajorityRule)


# Sample size
Data_ox <- dataGenus %>% 
  mutate(Oxygen_trait = Oxygen_Conservative) %>%  
  dplyr::select(Oxygen_trait,zscore,Genus,Genus_Phylo,Phylum,Phylum_simplified) %>% 
  drop_na()
dim(Data_ox)[1]

# PGLS Model
nTree=length(myGenusTrees)
Pgls_outupt = tibble(F_val=rep(NA,nTree),p=rep(NA,nTree),lambda=rep(NA,nTree))

for (i in 1:nTree) {
  dataCaper<- comparative.data(myGenusTrees[[i]],  Data_ox, Genus, vcv=TRUE, vcv.dim=3)
  PGLS <- pgls(zscore ~ Oxygen_trait, data = dataCaper, lambda='ML')
  aov_PGLS =  anova.pgls(PGLS) 
  
  Pgls_outupt[i,"F_val"] = aov_PGLS$`F value`[1]
  Pgls_outupt[i,"p"] = aov_PGLS$`Pr(>F)`[1]
  Pgls_outupt[i,"lambda"] = PGLS$param["lambda"]
}

median_PGLS_output = apply(Pgls_outupt,2,median)
median_PGLS_output 
sd_PGLS_output = apply(Pgls_outupt,2,sd)
sd_PGLS_output


if (median_PGLS_output["p"]<.001) {kp="<.001"} else {kp=round(median_PGLS_output["p"],3)}


PGLS_model_text = paste("F value=",round(median_PGLS_output["F_val"],2),
                        "\np=",kp,
                        "; n=",dim( Data_ox)[1],
                        sep="")


OxygenPlotY = Data_ox %>% 
  subset(!is.na(Oxygen_trait)) %>% 
  ggplot(aes(y=-zscore,
             x=Oxygen_trait,na.rm = T))+ 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(col=Phylum_simplified),alpha=alphaJitter,size=sizeJitter) + scale_color_manual(values=cols)  +
  ylab(TaxoSpecificitynameY) + xlab("Oxygen requirement")+
  annotate(geom="text", x=1.5,y=Ystats , label=PGLS_model_text,size=P2*0.35) +
  ylim(c(YlimDown,YlimUp )) +
  MyTheme + theme(legend.position="none")+ ggtitle("YEA dataset")

OxygenPlotY

# Combine plots
## ----------

fig6=plot_grid(OxygenPlotS,OxygenPlotY,
               nrow = 1, 
               align="h",axis = "l")
fig6
ggsave(plot = fig6,filename = outputFig ,height = 70,width = 180,units = "mm")




