##########################################
# Supp Fig. Bacterial TREE    
##########################################
library(ggtree)

# Run the Prep fig3-4 scripts first. 
# Plot one example of Bacterial Phylogeny 

output_file = "Submission/Mol_Ecol/V2/Supp_data/Supp_Figure_5.pdf"

tree = myGenusTrees[[25]]

dataGenus = Genus_Specificity %>% 
  rename(zscore=sr.obs.z_SONG_morethan2_samples_per_host_5k_median,
         zscore_altforPD=SR_Z_PDsampling_median,
         zscorePD=sesPD_withoutZero_median) %>% 
  dplyr::select(Genus,Genus_Phylo,n, Phylum_simplified,Phylum, zscore,zscorePD,zscore_altforPD,BSD,
                Oxygen_Conservative,Oxygen_MajorityRule,Spore_Conservative,Spore_MajorityRule,ASVs_in_Genus_SONG_morethan2_samples_per_host_5k_sum)


# Plots

p <- ggtree(tree) + 
      xlim(0, 1.2) + # to allow more space for labels
  geom_treescale(x=.8) 


p <- p %<+% dplyr::select(dataGenus,Genus,Phylum) + 
  geom_tiplab(aes(fill = factor(Phylum)),
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.15, "lines"), # amount of padding around the labels
              label.size = 0,
              size=1.5)  +# size of label border
  theme(legend.position = c(0.8,0.8), 
        legend.title = element_blank(), # no title
        legend.key = element_blank()) # no keys


ggsave(filename = output_file,plot = p,,height = 220,width = 120,units = "mm")


