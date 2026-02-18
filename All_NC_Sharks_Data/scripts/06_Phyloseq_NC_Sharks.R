#use Phyloseq to visualize NC_Sharks_Data
#Author Eldridge Wisely

#Load Libraries

library(ggplot2)

#Load phyloseq object of just fecal swabs made in 05 script:
NC_Sharks.ps<-readRDS("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/NC_Sharks.ps.RDS")
#Starting Phlyoseq!

NC_Sharks.ps
sample_names(NC_Sharks.ps)
rank_names(NC_Sharks.ps)
sample_variables(NC_Sharks.ps)

NC_Sharks.ps

plot_richness(NC_Sharks.ps, x="material","Species..scientific.name.", measures=c("Shannon", "Simpson"), color="material")

plot_richness(NC_Sharks.ps, x="Species..scientific.name.","material", measures=c("Shannon", "Simpson"), color="Species..scientific.name.")

ggplot2::ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/alpha_diversity_measures_by_species_and_material.jpg")

# Transform data to proportions as appropriate for Bray-Curtis distances
NC_Sharks.ps.prop <- transform_sample_counts(NC_Sharks.ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(NC_Sharks.ps.prop, method="NMDS", distance="bray")

#Warning message:
#In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data

plot_ordination(NC_Sharks.ps.prop, ord.nmds.bray, color="Species..scientific.name."
                , title="Bray NMDS")
#changed color to species.

top20 <- names(sort(taxa_sums(NC_Sharks.ps), decreasing=TRUE))[0:20]
ps.top20 <- transform_sample_counts(NC_Sharks.ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Species..scientific.name.",fill="genus")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")
ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/Bar_Chart_of_Species_by_top20preyGenus.jpg")
###############

#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(NC_Sharks.ps))
standf = function(x, t=total) round(t * (x / sum(x)))
NC_Sharks.ps.med.norm = transform_sample_counts(NC_Sharks.ps, standf)

plot_bar(NC_Sharks.ps.med.norm, x="Species..scientific.name.",fill="class")+ 
  geom_bar(aes(color=class, fill=class), stat="identity", position="stack")

ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/Bar_Chart_of_Species_by_PreyClass.jpg")

plot_bar(NC_Sharks.ps.med.norm, fill = "genus", facet_grid = "Species..scientific.name.")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")
ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/Faceted_sample_species_prey-genus.jpg")

plot_bar(NC_Sharks.ps.med.norm, fill = "family", facet_grid = "Species..scientific.name.")+ 
  geom_bar(aes(color=family, fill=family), stat="identity", position="stack")
ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/Faceted_sample_species_prey-family.jpg")

#Merge samples by donor species
NC_Shark_donor_merge.ps <- merge_samples(NC_Sharks.ps, "Species..scientific.name.")
plot_bar(NC_Shark_donor_merge.ps, fill = "family") + 
  geom_bar(aes(color=family, fill=family), stat="identity", position="stack")

#different organization of the graph
plot_bar(NC_Sharks.ps.med.norm, x="genus", fill = "genus", facet_grid = Species..scientific.name.~material) +
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

#ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/Stomachs_vs_Fecals_for_Each_Species_by_prey_genus.jpg")

#just the top 20
plot_bar(ps.top20, x="genus", fill = "genus", facet_grid = Species..scientific.name.~material) +
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")


NC_Sharks.ps.ord <- ordinate(NC_Sharks.ps.prop, "NMDS", "bray")
plot_ordination(NC_Sharks.ps,NC_Sharks.ps.ord , type="material", color="Species..scientific.name.", 
                shape="material", title="biplot", label = "material") +  
  geom_point(size=3)



#### Make a Final Taxa Table with the number of reads assigned to each genus or species ####

#This is the most important part to me to make the otu table readable.
#merge motus by prey species
NC_Shark_prey_species_merge.ps <- tax_glom(NC_Sharks.ps, taxrank = "species", NArm = FALSE)
plot_bar(NC_Shark_prey_species_merge.ps, fill = "family") + 
  geom_bar(aes(color=family, fill=family), stat="identity", position="stack")



NC_Shark_prey_species_merge.ps
sample_variables(NC_Shark_prey_species_merge.ps)
rank_names(NC_Shark_prey_species_merge.ps)

#export read table with species or genus (if only genus available) instead of motu/asv now that they've been properly merged
taxtbl <- as.data.frame(tax_table(NC_Shark_prey_species_merge.ps))
otutbl<-as.data.frame(otu_table(NC_Shark_prey_species_merge.ps))
sampletbl<-as.data.frame(sample_data(NC_Shark_prey_species_merge.ps))

write.csv(otutbl, "All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/preyspeciesmergeotutable.csv")
#looks good
library(tibble)
library(dplyr)

otutbl$species <-
  taxtbl$species[match(rownames(otutbl),rownames(taxtbl))]

otutbl$genus <-
  taxtbl$genus[match(rownames(otutbl),rownames(taxtbl))]

df<-otutbl %>% relocate(genus)
df<-df %>% relocate(species, .after = genus)

df1<- as.data.frame(t(df))

df1<- rownames_to_column(df1, var="sample")
sampletbl.1<-rownames_to_column(sampletbl, var="sample")
df2<- left_join(df1, sampletbl.1, by ="sample")

colnames(df2)

df3<- df2 %>% select(-"type",-"plate_no",-"plate_col.x",-"plate_row",-"tag_fwd",-"tag_rev",-"primer_fwd.x",-"primer_rev.x",-"project",-"control_type",-"nb_reads.x",-"nb_motus.x",-"low_contamination_level.x",-"seqdepth_ok.x",-"nb_reads_postmetabaR.x",-"nb_motus_postmetabaR.x",-"plate_col.y",-"plate_col.y",-"primer_fwd.y",-"primer_rev.y",-"nb_reads.y",-"nb_motus.y",-"low_contamination_level.y",-"seqdepth_ok.y",-"nb_reads_postmetabaR.y",-"nb_motus_postmetabaR.y")
colnames(df3)

write.csv(df3, "All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/NC_Sharks_samples_by_species_combined_table.csv")

#### Done making final taxa table ####
library(microViz)

NC_Shark_prey_species_merge.ps<-NC_Shark_prey_species_merge.ps %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )

saveRDS(NC_Shark_prey_species_merge.ps, file= "All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/taxa_named_NC_Shark_prey_species_merge.ps.RDS")


plot_richness(NC_Shark_prey_species_merge.ps, x="Species..scientific.name.","material", measures=c("Shannon", "Simpson"), color="Species..scientific.name.")

ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/alpha_diversity_measures_by_species_and_material.jpg")

top20species <- names(sort(taxa_sums(NC_Shark_prey_species_merge.ps), decreasing=TRUE))[0:20]
ps.top20species <- transform_sample_counts(NC_Shark_prey_species_merge.ps, function(OTU) OTU/sum(OTU))
ps.top20species <- prune_taxa(top20species, ps.top20species)
plot_bar(ps.top20species, x="Species..scientific.name.",fill="genus")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/Bar_Chart_of_Species_by_top20preyGenus.jpg")


plot_bar(ps.top20species, x="Species..scientific.name.",fill="species", )+ 
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")

ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/Bar_Chart_of_Species_by_top20preySpecies_incl_NA.jpg")

#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(NC_Shark_prey_species_merge.ps))
standf = function(x, t=total) round(t * (x / sum(x)))
NC_Sharks.ps.species.med.norm = transform_sample_counts(NC_Shark_prey_species_merge.ps, standf)

#different organization of the graph
plot_bar(NC_Sharks.ps.species.med.norm, x="genus", fill = "genus", facet_grid = Species..scientific.name.~material) +
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/Stomachs_vs_Fecals_for_Each_Species_by_prey_genus1.jpg")

plot_bar(NC_Sharks.ps.species.med.norm, fill = "genus", facet_grid = "Species..scientific.name.")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")
ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/Faceted_sample_species_prey-genus.jpg")

plot_bar(NC_Sharks.ps.species.med.norm, x="Species..scientific.name.",fill="class")+ 
  geom_bar(aes(color=class, fill=class), stat="identity", position="stack")

ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/Bar_Chart_of_Species_by_PreyClass.jpg")

NC_Sharks.ps.species.med.norm






#the figure you wanted..
NC_Shark_prey_species_merge.ps.prop <- transform_sample_counts(NC_Shark_prey_species_merge.ps, function(species) species/sum(species))

plot_bar(NC_Shark_prey_species_merge.ps.prop,x="sample_id", fill="genus") + 
  facet_grid(~Species..scientific.name., scales="free") +
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")
ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/read_proportions_by_shark_species_prey_genus.jpg")


plot_heatmap(NC_Shark_prey_species_merge.ps.prop, method = "NMDS", distance = "bray", taxa.label = "species", sample.label="sample_id")
ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/heatmap_sharkID_by_species_proportions.jpg")


plot_bar(NC_Shark_prey_species_merge.ps.prop, x="genus", fill = "genus", facet_grid = Species..scientific.name.~material) +
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

plot_bar(NC_Shark_prey_species_merge.ps.prop, fill="genus") + facet_wrap(~Species..scientific.name., scales= "free_x", nrow=1) +
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/read_proportions_by_shark_species_prey_genus.jpg")

plot_heatmap(NC_Shark_prey_species_merge.ps.prop, method = "NMDS", distance = "bray", taxa.label = "genus" )
ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/heatmap_sample_by_prey_genus_proportions.jpg")

plot_richness(NC_Shark_prey_species_merge.ps, x="Fork.Length..FL..mm..", color="Species..scientific.name.", measures=c("Observed"))

ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/alpha_diversity_prey_richness_by_fork_length_all_species.jpg")

plot_richness(NC_Shark_prey_species_merge.ps, x="Fork.Length..FL..mm..", color="Species..scientific.name.", measures=c("Observed", "Shannon")) + geom_boxplot()

plot_richness(NC_Shark_prey_species_merge.ps, x="Species..scientific.name.", color="Species..scientific.name.", measures=c("Observed", "Shannon")) + geom_boxplot()
ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/prey_species_diversity_obs_Shannon_by_Shark_species.jpg")

#following https://alexiscarter.github.io/metab/Phyloseq_script.html

sample_data(NC_Sharks.ps)$Species = factor(sample_data(NC_Sharks.ps)$Species..scientific.name., levels = c('Rhizoprionodon terraenovae', 'Carcharhinus limbatus', 'Carcharhinus acronotus','Sphyrna tiburo'))
#sample_data(NC_Sharks.ps)$type = factor(sample_data(NC_Sharks.ps)$material, levels = c('fecal', 'stomach'))

NC_Sharks.ps.log <- transform_sample_counts(NC_Sharks.ps, function(x) log(x+1))

NC_Sharks.ps.log.sp = tax_glom(NC_Sharks.ps.log, "species", NArm = FALSE)

ntaxa(NC_Sharks.ps.log)
nsamples(NC_Sharks.ps.log)
rank_names(NC_Sharks.ps.log)
sample_variables(NC_Sharks.ps.log)
otu_table(NC_Sharks.ps.log)
tax_table(NC_Sharks.ps.log)

plot_bar(NC_Sharks.ps.log.sp, fill = "family") +
  geom_bar(stat="identity") +
  facet_wrap(~Species, ncol = 4, scales = "free_x") +
  labs(x = "sample", y = "Reads abundance of prey")

plot_bar(NC_Sharks.ps.log.sp, fill = "genus") +
  geom_bar(stat="identity") +
  facet_wrap(type~Species, ncol = 4, scales = "free_x") +
  labs(x = "sample", y = "log normalized read abundance")


#Following the microviz tutorial https://david-barnett.github.io/microViz/index.html
#awesome stuff!

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)

#install.packages(
#  "microViz",
#  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
#)

#install.packages("ggraph") # for taxatree_plots()
#install.packages("DT") # for tax_fix_interactive()
#install.packages("corncob") # for example datasets and beta binomial models
library(microViz)
library(ggraph)
library(DT)
library(corncob)

NC_Shark_prey_taxa_merge.ps<-tax_glom(NC_Sharks.ps, taxrank = "species",NArm=FALSE)

#tax_fix_interactive(NC_Shark_prey_taxa_merge.ps)
NC_Shark_prey_taxa_merge.ps<-NC_Shark_prey_taxa_merge.ps %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )

top20species1 <- names(sort(taxa_sums(NC_Shark_prey_taxa_merge.ps), decreasing=TRUE))[0:20]
ps.top20species1 <- transform_sample_counts(NC_Shark_prey_taxa_merge.ps, function(OTU) OTU/sum(OTU))
ps.top20species1 <- prune_taxa(top20species1, ps.top20species1)
plot_bar(ps.top20species1, x="Species..scientific.name.",fill="species")+ 
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")

ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/stacked_bar_chart_top20species_level_prey_by_shark_species.jpg")

NC_Shark_prey_taxa_merge.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 34 taxa and 80 samples ]
# sample_data() Sample Data:       [ 80 samples by 62 sample variables ]
# tax_table()   Taxonomy Table:    [ 34 taxa by 7 taxonomic ranks ]

top34species1 <- names(sort(taxa_sums(NC_Shark_prey_taxa_merge.ps), decreasing=TRUE))[0:34]
ps.top34species1 <- transform_sample_counts(NC_Shark_prey_taxa_merge.ps, function(OTU) OTU/sum(OTU))
ps.top20species1 <- prune_taxa(top43species1, ps.top43species1)
plot_bar(ps.top34species1, x="Species..scientific.name.",fill="species")+ 
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")

#almost too many species for the legend to make sense

# play with ordinations interactively
pseq <- NC_Shark_prey_taxa_merge.ps %>%
  phyloseq_validate()
#uncomment the next line for a web browser based way to explore ordinations with these data!
#ord_explore(pseq) 
#insert code from the above program that you want to keep... 

#PCA of binary transformed (presence/absence) prey genus.
pseq %>%
  tax_transform(rank = "genus", trans = "binary") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Species", fill = "Species",
    shape = "circle", alpha = 0.5,
    size = 2
  ) +
  ggplot2::stat_ellipse(
    ggplot2::aes(colour = Species)
  )

ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/presence_absence_diet_PCA.jpg")

#compositional transformation of data, showing arrows of the differentiating genuses and the strength of their influence on the PCA plot
pseq %>%
  tax_transform(rank = "species", trans = "compositional") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    plot_taxa = 1:3,
    colour = "Species..scientific.name.", fill = "Species..scientific.name.",
    shape = "circle", alpha = 0.5,
    size = 2
  )

ggsave("All_NC_Sharks_Data/data/07_Phyloseq_results_and_Visualizations/PCA_of_diet_composition_with_arrows.jpg")


