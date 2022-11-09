
myGenusTrees = readRDS("Data/MicrobialGenusTrees.RDS")

##########################################################
######        Results  for PHYLO specificity        ######
##########################################################


datasets <- c("SEA","YEA")
Explanatory_var <- c("PD","BSD","Oxygen","Spore")

Stats <- expand.grid(dataset=datasets,var=Explanatory_var,
                     F_factor=NA,pval_model=NA,
                     Lambda_modelResiduals=NA,
                     stringsAsFactors = F)


for (data in c("SEA","YEA")) {
  
  if (data=="SEA") {
    ## 2. Supp plot: Song et phylo range
    ## ----------------------------------## 
    outputFig1 <- "Submission/Mol_Ecol/V2/Supp_data/Supp_Figure_2.pdf"
    dataSp=ASV_specificity_Song %>% 
      rename(zscore=sesPD_SONG_morethan2_samples_per_host_5k)
    
    outputFig2 <- "Submission/Mol_Ecol/V2/Supp_data/Supp_Figure_9.pdf"
    dataGenus = Genus_Specificity %>% 
      rename(zscore=sesPD_withoutZero_median) %>% 
      dplyr::select(Genus,Genus_Phylo,n, Phylum_simplified,Phylum, zscore,BSD,
             Oxygen_Conservative,Oxygen_MajorityRule,Spore_Conservative,Spore_MajorityRule,ASVs_in_Genus_SONG_morethan2_samples_per_host_5k_sum)
    
    
  }
  
  if (data=="YEA"){
    ## 4. Supp plot: Younglbut et phylo range
    outputFig1 <- "Submission/Mol_Ecol/V2/Supp_data/Supp_Figure_4.pdf"
    dataSp=ASV_specificity_Y %>% 
      rename(zscore=sesPD_morethan2_samples_per_host_5k_Youngblut)
    
    outputFig2 <- "Submission/Mol_Ecol/V2/Supp_data/Supp_Figure_10.pdf"
    dataGenus = Genus_Specificity %>% 
      rename(zscore=sesPD_morethan2_samples_per_host_5k_Youngblut_median) %>% 
      dplyr::select(Genus,Genus_Phylo,n, Phylum_simplified,Phylum, zscore,BSD,
             Oxygen_Conservative,Oxygen_MajorityRule,Spore_Conservative,Spore_MajorityRule)
    
    
  }  
  



  

# --------------------#
# -----------------
# AXis name 
# -----------------
# --------------------#

TaxoSpecificityname = "Taxonomic Specificity (Z-score)"
PhyloSpecificityname  = "Phylogenetic Specificity (Z-score)"



# --------------------#
# -----------------
# ASV level figure 
# -----------------
# --------------------#

# Histograms on specificity 
## ------------------------
Specificity_hist_plot = dataSp %>% 
  ggplot(aes(y=-zscore))+ 
  geom_histogram() + xlab("")+
  MyTheme + ylab(PhyloSpecificityname)+geom_hline(yintercept=1.96, linetype="dashed") + ggtitle("A")
Specificity_hist_plot 


# Genus plot 
## ----------------

# Sample size
s <- dataSp %>% 
  dplyr::select(zscore,Genus) %>% 
  drop_na() %>% 
  group_by(Genus) %>% 
  summarise(n=n()) 

# Test
kruskalResult = kruskal.test(zscore~Genus,dataSp)
kruskalResult
if (kruskalResult$p.value<.001) {kp="<.001"} else {kp=round(kruskalResult$p.value,3)}
kruskatest=paste("KW statistic=",round(kruskalResult$statistic,1),
                 "; p value=",kp,
                 "; n(ASVs)=",sum(s$n),
                 "; n(genus)=",dim(s)[1],
                 sep="")
kruskatest

# Plot
Genus_aggregated_plot = dataSp %>% 
  ggplot(aes(y=-zscore,
             x=reorder(Genus,zscore,
                       median,na.rm = TRUE)))+ 
  geom_boxplot(position=position_dodge(width=13.8), outlier.size = .3,outlier.colour = "grey",aes(col=Phylum_simplified))+geom_hline(yintercept=1.96, linetype="dashed") +
  ylab("") +xlab("Genus") +
  scale_color_manual(values=cols) +
  annotate(geom="text", x=30,y=-2, label=kruskatest,size=P2*0.35) +
  MyTheme + ggtitle("B") + theme(legend.title=element_blank(), axis.text.x =element_text(size=P1))

legend <- cowplot::get_legend(Genus_aggregated_plot)

Genus_aggregated_plot  <- Genus_aggregated_plot  + theme(legend.position="none")
Genus_aggregated_plot 

# Combine plots
## ------------

fig5=plot_grid(Specificity_hist_plot ,Genus_aggregated_plot ,
               nrow = 1, 
               rel_widths = c(.2,1),
               align="h",axis = "l")

fig5=plot_grid(fig5 ,legend, nrow = 1, 
               rel_widths = c(1,.07))

ggsave(plot = fig5,filename = outputFig1 ,height = 120,width = 360,units = "mm")



# --------------------#
# --------------------#
# Genus level figure 
# --------------------#
# --------------------#

Ystats = 13.5
YlimUp = 14
YlimDown=-1
alphaJitter = 1
sizeJitter =.5

## ------------- ##
# LINK TO BSD  
## ------------  ##

Variable="BSD"
outliers=c("Fusobacterium","Ureaplasma")

# Sample size
BSD_data <- dataGenus %>% 
  dplyr::select(BSD,zscore,Genus,Genus_Phylo,Phylum_simplified) %>% 
  subset(!Genus%in%outliers) %>% 
  drop_na()

# PGLS Model
nTree=length(myGenusTrees)
Pgls_outupt = tibble(F_val=rep(NA,nTree),p=rep(NA,nTree),lambda=rep(NA,nTree))

for (i in 1:nTree) {
  dataCaper<- comparative.data(myGenusTrees[[i]],  BSD_data, Genus, vcv=TRUE, vcv.dim=3)
  PGLS <- try(pgls(BSD ~ zscore, data = dataCaper, lambda='ML'))
  if (!class(PGLS)=="try-error") {
    aov_PGLS =  anova.pgls(PGLS) 
    
    Pgls_outupt[i,"F_val"] = aov_PGLS$`F value`[1]
    Pgls_outupt[i,"p"] = aov_PGLS$`Pr(>F)`[1]
    Pgls_outupt[i,"lambda"] = PGLS$param["lambda"]
  }

}

median_PGLS_output = apply(Pgls_outupt,2,median,na.rm=T)
median_PGLS_output 
sd_PGLS_output = apply(Pgls_outupt,2,sd,na.rm=T)
sd_PGLS_output

Stats <- Stats %>% mutate(Lambda_modelResiduals=ifelse(dataset==data&var==Variable,  median_PGLS_output["lambda"],Lambda_modelResiduals),
                          F_factor=ifelse(dataset==data&var==Variable,median_PGLS_output["F_val"],F_factor),
                          pval_model=ifelse(dataset==data&var==Variable,median_PGLS_output["p"],pval_model))

PGLS_model_text = paste("F value=",round(median_PGLS_output["F_val"],2),
                        "\np=",round(median_PGLS_output["p"],3),
                        "; n=",dim(BSD_data)[1],
                        sep="")

# PLOT
BSDplotGenus_Specificity = BSD_data %>% 
  ggplot(aes(y=-zscore,x=BSD))+   
  geom_point(aes(col=Phylum_simplified)) + scale_color_manual(values=cols) + theme(legend.position="none") +
  geom_smooth(col="black",se = F,method = "lm")+
  ylab(PhyloSpecificityname) + xlab("Transmission-mode score")+
  #annotate(geom="text", x=1.2,y=Ystats , label=PGLS_model_text ,size=P2*0.35) +
  ylim(c(YlimDown,YlimUp )) +
  MyTheme + theme(legend.position="none") + ggtitle("B")

BSDplotGenus_Specificity





## ------------- ##
# LINK TO OXYGEN 
## -------------

Variable="Oxygen"

# Data
Data_ox <- dataGenus %>% 
  #mutate(Oxygen_trait = Oxygen_Conservative) %>%  
  mutate(Oxygen_trait = Oxygen_MajorityRule) %>%  
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


Stats <- Stats %>% mutate(Lambda_modelResiduals=ifelse(dataset==data&var==Variable,  median_PGLS_output["lambda"],Lambda_modelResiduals),
                          F_factor=ifelse(dataset==data&var==Variable,median_PGLS_output["F_val"],F_factor),
                          pval_model=ifelse(dataset==data&var==Variable,median_PGLS_output["p"],pval_model))

if (median_PGLS_output["p"]<.001) {kp="<.001"} else {kp=round(median_PGLS_output["p"],3)}


PGLS_model_text = paste("F value=",round(median_PGLS_output["F_val"],2),
                        "\np=",kp,
                        "; n=",dim( Data_ox)[1],
                        sep="")


OxygenPlot = Data_ox %>% 
  subset(!is.na(Oxygen_trait)) %>% 
  ggplot(aes(y=-zscore,
             x=Oxygen_trait,na.rm = T))+ 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(col=Phylum_simplified),alpha=alphaJitter,size=sizeJitter) + scale_color_manual(values=cols)  +
  ylab(PhyloSpecificityname) + xlab("Oxygen requirement")+
  annotate(geom="text", x=1.5,y=Ystats , label= PGLS_model_text,size=P2*0.35) +
  ylim(c(YlimDown,YlimUp )) +
  MyTheme + theme(legend.position="none")+ ggtitle("C")

OxygenPlot 

## ------------- ##
# LINK TO SPORE
## ------------- ##

Variable="Spore"
# Sample size
dataSpore <- dataGenus %>% 
  subset(Phylum=="Firmicutes") %>% 
  mutate(spore_trait = Spore_Conservative) %>% 
  dplyr::select(spore_trait,zscore,Genus,Genus_Phylo,Phylum_simplified) %>% 
  drop_na()
dim(dataSpore)[1]

# PGLS Model
nTree=length(myGenusTrees)
Pgls_outupt = tibble(F_val=rep(NA,nTree),p=rep(NA,nTree),lambda=rep(NA,nTree))

for (i in 1:nTree) {
  dataCaper<- comparative.data(myGenusTrees[[i]],  dataSpore, Genus, vcv=TRUE, vcv.dim=3)
  PGLS <- pgls(zscore ~ spore_trait, data = dataCaper, lambda='ML')
  aov_PGLS =  anova.pgls(PGLS) 
  
  Pgls_outupt[i,"F_val"] = aov_PGLS$`F value`[1]
  Pgls_outupt[i,"p"] = aov_PGLS$`Pr(>F)`[1]
  Pgls_outupt[i,"lambda"] = PGLS$param["lambda"]
}

median_PGLS_output = apply(Pgls_outupt,2,median)
median_PGLS_output 
sd_PGLS_output = apply(Pgls_outupt,2,sd)
sd_PGLS_output


Stats <- Stats %>% mutate(Lambda_modelResiduals=ifelse(dataset==data&var==Variable,  median_PGLS_output["lambda"],Lambda_modelResiduals),
                          F_factor=ifelse(dataset==data&var==Variable,median_PGLS_output["F_val"],F_factor),
                          pval_model=ifelse(dataset==data&var==Variable,median_PGLS_output["p"],pval_model))

PGLS_model_text = paste("F value=",round(median_PGLS_output["F_val"],2),
                        "\np=",round(median_PGLS_output["p"],3),
                        "; n=",dim(dataSpore)[1],
                        sep="")


SporePlot = dataSpore %>% 
  subset(!is.na(spore_trait)) %>% 
  ggplot(aes(y=-zscore,
             x=spore_trait,na.rm = T))+ 
  geom_boxplot() +
  geom_jitter(aes(col=Phylum_simplified),alpha=alphaJitter,size=sizeJitter) + scale_color_manual(values=cols)  +
  ylab(PhyloSpecificityname) + xlab("Spore formation abilities")+
  annotate(geom="text", x=1.5,y=Ystats, label=PGLS_model_text,size=P2*0.35) +
  ylim(c(YlimDown,YlimUp )) +
  MyTheme + theme(legend.position="none") + ggtitle("D")

SporePlot


# Combine plots
## ----------

fig6=plot_grid(BSDplotGenus_Specificity,OxygenPlot,SporePlot,
               nrow = 1, 
               align="h",axis = "l")
fig6=plot_grid(fig6 ,legend, nrow = 1, 
               rel_widths = c(1,.2))
fig6
ggsave(plot = fig6,filename = outputFig2 ,height = 70,width = 180,units = "mm")




}

Stats 

