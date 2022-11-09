##########################################
# Supp Fig. Comparisons across datasets  
##########################################


# Comparaison of microbial specificity estimtes across datasets
# Run the Prep fig3-4 scripts first. 


# plot params 
Ystats = 13.5
YlimUp = 14
YlimDown=-1
alphaJitter = 1
sizeJitter =.5

outputFig <- "Submission/Mol_Ecol/V2/Supp_data/Supp_Figure_6.pdf"



#data
dataGenusCompa = Genus_Specificity %>% 
  rename(zscoreSong=sr.obs.z_SONG_morethan2_samples_per_host_5k_median,
         zscorePDSong=sesPD_withoutZero_median,
         zscoreY=sr.obs.z_morethan2_samples_per_host_5k_Youngblut_median,
         zscorePDY=sesPD_morethan2_samples_per_host_5k_Youngblut_median) %>% 
  dplyr::select(Genus,n,Genus_Phylo, Phylum_simplified,Phylum, zscoreSong,zscorePDSong,zscoreY,zscorePDY) %>% 
  drop_na()

# # # # # # # # 
# host range ##
# # # # # # # # 


# PGLS Model
nTree=length(myGenusTrees)
Pgls_outupt = tibble(F_val=rep(NA,nTree),p=rep(NA,nTree),lambda=rep(NA,nTree))

for (i in 1:nTree) {
  dataCaper<- comparative.data(myGenusTrees[[i]], dataGenusCompa, Genus, vcv=TRUE, vcv.dim=3)
  PGLS <- pgls(zscoreSong ~ zscoreY, data = dataCaper, lambda='ML')
  aov_PGLS =  anova.pgls(PGLS) 
  
  Pgls_outupt[i,"F_val"] = aov_PGLS$`F value`[1]
  Pgls_outupt[i,"p"] = aov_PGLS$`Pr(>F)`[1]
  Pgls_outupt[i,"lambda"] = PGLS$param["lambda"]
}

median_PGLS_output = apply(Pgls_outupt,2,median)
median_PGLS_output
sd_PGLS_output = apply(Pgls_outupt,2,sd)
sd_PGLS_output

PGLS_model_text = paste("F value=",round(median_PGLS_output["F_val"],2),
      "\np=",round(median_PGLS_output["p"],3),
      "; n=",dim(dataGenusCompa )[1],
      sep="")


CompaPlotSR = dataGenusCompa %>% 
  ggplot(aes(y=-zscoreSong,
             x=-zscoreY,na.rm = T))+ 
  geom_smooth(col="black",se = F,method = "lm")+
  ylab("Taxo. Specificity (Z-score) SEA dataset") + xlab("taxo. Specificity (Z-score) YEA dataset")+
  MyTheme +   annotate(geom="text", x=1.5,y=Ystats , label=PGLS_model_text,size=P2*0.35)+
  geom_point(aes(col=Phylum_simplified)) + scale_color_manual(values=cols) 


CompaPlotSR 
ggsave(plot = CompaPlotSR ,filename = outputFig ,height = 120,width = 120,units = "mm")





