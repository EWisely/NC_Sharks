#use Phyloseq to visualize NC_Sharks_Data
#Author Eldridge Wisely

library(ggplot2)
library(phyloseq)
library (vegan)
library(metagMisc)

#load phyloseq object
NC_Sharks.ps<-readRDS(file="data/Paired/07_Phyloseq/Paired_NC_Sharks.ps")
NC_Sharks.ps

sample_names(NC_Sharks.ps)
rank_names(NC_Sharks.ps)
sample_variables(NC_Sharks.ps)


plot_richness(NC_Sharks.ps, x="material","Species..scientific.name.", measures=c("Shannon", "Simpson"), color="material")

plot_richness(NC_Sharks.ps, x="Species..scientific.name.","material", measures=c("Shannon", "Simpson"), color="Species..scientific.name.")

ggplot2::ggsave("data/Paired/07_Phyloseq/alpha_diversity_measures_by_species_and_material.jpg")

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
ggsave("data/Paired/07_Phyloseq/Bar_Chart_of_Species_by_top20preyGenus.jpg")
###############

#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(NC_Sharks.ps))
standf = function(x, t=total) round(t * (x / sum(x)))
NC_Sharks.ps.med.norm = transform_sample_counts(NC_Sharks.ps, standf)

plot_bar(NC_Sharks.ps.med.norm, x="Species..scientific.name.",fill="class")+ 
  geom_bar(aes(color=class, fill=class), stat="identity", position="stack")

ggsave("data/Paired/07_Phyloseq/Bar_Chart_of_Species_by_PreyClass.jpg")

plot_bar(NC_Sharks.ps.med.norm, fill = "genus", facet_grid = "Species..scientific.name.")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")
ggsave("data/Paired/07_Phyloseq/Faceted_sample_species_prey-genus.jpg")

plot_bar(NC_Sharks.ps.med.norm, fill = "family", facet_grid = "Species..scientific.name.")+ 
  geom_bar(aes(color=family, fill=family), stat="identity", position="stack")
ggsave("data/Paired/07_Phyloseq/Faceted_sample_species_prey-family.jpg")


#Merge samples by donor species
NC_Shark_donor_merge.ps <- merge_samples(NC_Sharks.ps, "Species..scientific.name.")
plot_bar(NC_Shark_donor_merge.ps, fill = "family") + 
  geom_bar(aes(color=family, fill=family), stat="identity", position="stack")

#different organization of the graph
plot_bar(NC_Sharks.ps.med.norm, x="genus", fill = "genus", facet_grid = Species..scientific.name.~material) +
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

ggsave("data/Paired/07_Phyloseq/Stomachs_vs_Fecals_for_Each_Species_by_prey_genus.jpg")

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

#transform to presence absence per species (and other taxonomic ranks as well)
NC_Shark_prey_species_merge_pa.ps<-phyloseq_standardize_otu_abundance(NC_Shark_prey_species_merge.ps, method = "pa")

#plot presence/absence data
plot_bar(NC_Shark_prey_species_merge_pa.ps, fill = "genus", x="Species..scientific.name.",facet_grid = "material" )+
  geom_bar(aes(color=genus, fill=genus), stat="identity", position = "stack")+
  ylab("Number of samples")
ggsave("data/Paired/07_Phyloseq/PA_paired_by_genus_stacked.jpg")

#now plot presence absence FOO by species without NAs
NC_Shark_prey_species_merge_sp.ps <- tax_glom(NC_Sharks.ps, taxrank = "species", NArm = TRUE)
NC_Shark_prey_species_merge_sp_pa.ps<-phyloseq_standardize_otu_abundance(NC_Shark_prey_species_merge_sp.ps, method = "pa")
plot_bar(NC_Shark_prey_species_merge_sp_pa.ps, fill = "species", x="Species..scientific.name.",facet_grid = "material" )+
  geom_bar(aes(color=species, fill=species), stat="identity", position = "stack")+
  ylab("Number of samples")
ggsave("data/Paired/07_Phyloseq/PA_paired_by_species_stacked.jpg")


#plot presence/absence data
plot_bar(NC_Shark_prey_species_merge_pa.ps, fill = "species", x="Species..scientific.name.",facet_grid = "material" )+
  geom_bar(aes(color=species, fill=species), stat="identity", position = "stack")+
  ylab("Number of samples")
ggsave("data/Paired/07_Phyloseq/PA_paired_by_species_stacked.jpg")


#different view of the presence/absence data
plot_bar(NC_Shark_prey_species_merge_pa.ps, x="genus", fill = "genus", facet_grid = Species..scientific.name.~material) +
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
  ylab("Number of samples")
ggsave("data/Paired/07_Phyloseq/PA_paired_by_genus.jpg")








#export read table with species or genus (if only genus available) instead of motu/asv now that they've been properly merged
taxtbl <- as.data.frame(tax_table(NC_Shark_prey_species_merge.ps))
otutbl<-as.data.frame(otu_table(NC_Shark_prey_species_merge.ps))
sampletbl<-as.data.frame(sample_data(NC_Shark_prey_species_merge.ps))


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

samples.df<-as.data.frame(NC_Shark_prey_species_merge.ps@sam_data)
samples.df.1<-rownames_to_column(samples.df, var="sample")
df2<- left_join(df1, samples.df.1, by ="sample")

colnames(df2)

df3<- df2 %>% select(-"type",-"plate_no",-"plate_col.x",-"plate_row",-"tag_fwd",-"tag_rev",-"primer_fwd.x",-"primer_rev.x",-"project",-"control_type",-"nb_reads.x",-"nb_motus.x",-"low_contamination_level.x",-"seqdepth_ok.x",-"nb_reads_postmetabaR.x",-"nb_motus_postmetabaR.x",-"date_collected.x",-"species.x",-"plate_col.y",-"primer_fwd.y",-"primer_rev.y",-"nb_reads.y",-"nb_motus.y",-"low_contamination_level.y",-"seqdepth_ok.y",-"nb_reads_postmetabaR.y",-"nb_motus_postmetabaR.y",-"date_collected.y",-"species.y")
colnames(df3)

write.csv(df3, "data/Paired/07_Phyloseq/NC_Sharks_paired_samples_by_species_combined_table.csv")

#### Done making final taxa table ####

#NC_Shark_prey_species_merge.ps
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 53 taxa and 131 samples ]
#sample_data() Sample Data:       [ 131 samples by 61 sample variables ]
#tax_table()   Taxonomy Table:    [ 53 taxa by 7 taxonomic ranks ]

#Make a subset with only paired fecals/stomachs


paired_ASVs<-as.data.frame(NC_Sharks.ps@tax_table)
library(tidyverse)

paired_ASVs<-paired_ASVs%>%
  mutate(taxon_name=
           coalesce(species,genus,family,order,class,phylum))

summary_paired_ASVs <- paired_ASVs |> 
  group_by(taxon_name) |> 
  summarize(number_of_ASVs = n())
mean(summary_paired_ASVs$number_of_ASVs)
median(summary_paired_ASVs$number_of_ASVs)
range(summary_paired_ASVs$number_of_ASVs)


write.csv(summary_paired_ASVs, file= "data/Paired/07_Phyloseq/summary_ASVs_per_taxon_paired.csv")


#NC_Shark_prey_species_merge_paired.ps<-subset_samples(NC_Shark_prey_species_merge.ps,Paired=="Yes")

NC_Shark_prey_species_merge_paired.ps<-NC_Shark_prey_species_merge.ps

plot_richness(NC_Shark_prey_species_merge.ps, x="Species..scientific.name.","material", measures=c("Shannon", "Simpson"), color="Species..scientific.name.")

ggsave("data/Paired/07_Phyloseq/alpha_diversity_measures_by_species_and_material.jpg")

top20species <- names(sort(taxa_sums(NC_Shark_prey_species_merge.ps), decreasing=TRUE))[0:20]
ps.top20species <- transform_sample_counts(NC_Shark_prey_species_merge.ps, function(OTU) OTU/sum(OTU))
ps.top20species <- prune_taxa(top20species, ps.top20species)
plot_bar(ps.top20species, x="Species..scientific.name.",fill="genus")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

ggsave("data/Paired/07_Phyloseq/Bar_Chart_of_Species_by_top20preyGenus.jpg")


plot_bar(ps.top20species, x="Species..scientific.name.",fill="species", )+ 
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")

ggsave("data/Paired/07_Phyloseq/Bar_Chart_of_Species_by_top20preySpecies_incl_NA.jpg")

#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(NC_Shark_prey_species_merge.ps))
standf = function(x, t=total) round(t * (x / sum(x)))
NC_Sharks.ps.species.med.norm = transform_sample_counts(NC_Shark_prey_species_merge.ps, standf)

#different organization of the graph
plot_bar(NC_Sharks.ps.species.med.norm, x="genus", fill = "genus", facet_grid = Species..scientific.name.~material) +
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

ggsave("data/Paired/07_Phyloseq/Stomachs_vs_Fecals_for_Each_Species_by_prey_genus1.jpg")

plot_bar(NC_Sharks.ps.species.med.norm, fill = "genus", facet_grid = "Species..scientific.name.")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")
ggsave("data/Paired/07_Phyloseq/Faceted_sample_species_prey-genus.jpg")

plot_bar(NC_Sharks.ps.species.med.norm, x="Species..scientific.name.",fill="class")+ 
  geom_bar(aes(color=class, fill=class), stat="identity", position="stack")

ggsave("data/Paired/07_Phyloseq/Bar_Chart_of_Species_by_PreyClass.jpg")

NC_Sharks.ps.species.med.norm



#Make Sorensen Similarity violin plot from paired fecal/stomach samples normalized to median sequencing depth... merged by prey species.

#NC_Shark_prey_species_merge_med.norm.paired.ps<-subset_samples(NC_Sharks.ps.species.med.norm,Paired=="Yes")
NC_Shark_prey_species_merge_med.norm.paired.ps<-NC_Sharks.ps.species.med.norm

NC_Shark_prey_species_merge_med.norm.paired.ps
#View(sample_data(NC_Shark_prey_species_merge_med.norm.paired.ps))

# Transform data to proportions as appropriate for Bray-Curtis distances
#NC_Sharks.ps.prop <- transform_sample_counts(NC_Sharks.ps, function(otu) otu/sum(otu))


#Eldridge's code for Sorensen similarity plots by species, genus, family and order

#Sørensen similarity plot to compare how similarity between the fecal and stomach samples from each individual
library(ggplot2)
library(phyloseq)


#View(NC_Shark_prey_species_merge_med.norm.paired.ps)
NC_Shark_prey_species_merge_med.norm.paired.ps
#View(sample_data(NC_Shark_prey_species_merge_med.norm.paired.ps))



#merge motus by prey species with no NAs
NC_Shark_prey_species_merge1.ps <- tax_glom(NC_Sharks.ps, taxrank = "species", NArm = TRUE)
plot_bar(NC_Shark_prey_species_merge1.ps, fill = "family") + 
  geom_bar(aes(color=family, fill=family), stat="identity", position="stack")

NC_Shark_prey_species_merge1.ps
#sample_variables(NC_Shark_prey_species_merge1.ps)

#after removing otus only ID'ed to genus level
#phyloseq-class experiment-level object (NArm=TRUE)
#otu_table()   OTU Table:         [ 43 taxa and 131 samples ]
#sample_data() Sample Data:       [ 131 samples by 61 sample variables ]
#tax_table()   Taxonomy Table:    [ 43 taxa by 7 taxonomic ranks ]


NC_Shark_prey_species_merge1.ps.prop <- transform_sample_counts(NC_Shark_prey_species_merge1.ps, function(species) species/sum(species))

#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(NC_Shark_prey_species_merge1.ps))
standf = function(x, t=total) round(t * (x / sum(x)))
NC_Sharks.ps.species1.med.norm = transform_sample_counts(NC_Shark_prey_species_merge1.ps, standf)


#subset the paired samples

#NC_Sharks.ps.species1.med.norm.paired.ps<-subset_samples(NC_Sharks.ps.species1.med.norm,Paired=="Yes")
NC_Sharks.ps.species1.med.norm.paired.ps<-NC_Sharks.ps.species1.med.norm

#extract matching fecal and stomach samples by Shark ID
library(tidyverse)



matched_fs_and_sts<- sample_data(NC_Sharks.ps.species1.med.norm.paired.ps)[,c("sample_id","material")]

matched_fs_and_sts$id<- rownames(sample_data(NC_Sharks.ps.species1.med.norm.paired.ps))
matched_fs_and_sts<-data.frame(matched_fs_and_sts)

matched_data_wide<-matched_fs_and_sts%>%
  pivot_wider(id_cols = sample_id,names_from = material, values_from = id)
View(matched_data_wide)

# Extract OTU table data for stomach and fecal samples
stomach_data <- otu_table(NC_Sharks.ps.species1.med.norm.paired.ps)[, matched_data_wide$stomach]
fecal_data <- otu_table(NC_Sharks.ps.species1.med.norm.paired.ps)[, matched_data_wide$fecal]

dim(stomach_data)
dim(fecal_data)

colnames(stomach_data)
colnames(fecal_data)

#calculate distances and put them in a symmetrical matrix
pairwise_sorensen_distances<-as.matrix(print(phyloseq::distance(NC_Sharks.ps.species1.med.norm.paired.ps, method="sor", binary=TRUE, type= "samples", upper=TRUE, diag=TRUE)))

View(pairwise_sorensen_distances)
class(pairwise_sorensen_distances)

#Yay, the matrix is symmetrical.  Now, I need to pull out the values per shark in matched_data_wide.

#make an empty column for numerical values
matched_data_wide$sor_dist<-NA_real_

for (shark_i in 1:nrow(matched_data_wide)){
  stomach_id<-matched_data_wide$stomach[shark_i]
  fecal_id<-matched_data_wide$fecal[shark_i]
  matched_data_wide$sor_dist[shark_i]<-pairwise_sorensen_distances[stomach_id,fecal_id]
}

# Plot similarity as a scatter plot
similarity_plot <- ggplot(data = matched_data_wide, aes(x = sample_id, y = sor_dist, fill = sor_dist)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +  # Adjust color scale
  labs(title = "Sørensen Pairwise Similarity between Prey Species in Stomach and Fecal Samples",
       x = "Shark_ID",
       y = "Sorensen Distances") +
  theme_minimal() +
  theme(text = element_text(size = 10),  # Adjust text size
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"))  # Adjust title size and style

print(similarity_plot)


#Do the same Sørensen distance matrix and scatterplot with genus tax_glom

#merge motus by prey genus
NC_Shark_prey_genus_merge.ps <- tax_glom(NC_Sharks.ps, taxrank = "genus", NArm = FALSE)
plot_bar(NC_Shark_prey_genus_merge.ps, fill = "genus") + 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

NC_Sharks.ps
NC_Shark_prey_species_merge.ps
NC_Shark_prey_genus_merge.ps


NC_Shark_prey_genus_merge.ps
sample_variables(NC_Shark_prey_genus_merge.ps)
rank_names(NC_Shark_prey_genus_merge.ps)
#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(NC_Shark_prey_genus_merge.ps))
standf = function(x, t=total) round(t * (x / sum(x)))
NC_Sharks.ps.genus.med.norm = transform_sample_counts(NC_Shark_prey_genus_merge.ps, standf)

NC_Sharks.ps.genus.med.norm

#Get just the fecal and stomach paired samples

#NC_Shark_prey_genus_merge_med.norm.paired.ps<-subset_samples(NC_Sharks.ps.genus.med.norm,Paired=="Yes")
NC_Shark_prey_genus_merge_med.norm.paired.ps<-NC_Sharks.ps.genus.med.norm
#re-do Sorensen distances for the genus phyloseq object

genus_matched_fs_and_sts<- sample_data(NC_Shark_prey_genus_merge_med.norm.paired.ps)[,c("sample_id","material")]

genus_matched_fs_and_sts$id<- rownames(sample_data(NC_Shark_prey_genus_merge_med.norm.paired.ps))
genus_matched_fs_and_sts<-data.frame(genus_matched_fs_and_sts)

genus_matched_data_wide<-genus_matched_fs_and_sts%>%
  pivot_wider(id_cols = sample_id,names_from = material, values_from = id)
#View(genus_matched_data_wide)

# Extract OTU table data for stomach and fecal samples
genus_stomach_data <- otu_table(NC_Shark_prey_genus_merge_med.norm.paired.ps)[, genus_matched_data_wide$stomach]
genus_fecal_data <- otu_table(NC_Shark_prey_genus_merge_med.norm.paired.ps)[, genus_matched_data_wide$fecal]

dim(genus_stomach_data)
dim(genus_fecal_data)

colnames(genus_stomach_data)
colnames(genus_fecal_data)

#calculate distances and put them in a symmetrical matrix
genus_pairwise_sorensen_distances<-as.matrix(print(phyloseq::distance(NC_Shark_prey_genus_merge_med.norm.paired.ps, method="sor", binary=TRUE, type= "samples", upper=TRUE, diag=TRUE)))

View(genus_pairwise_sorensen_distances)
class(genus_pairwise_sorensen_distances)

#Yay, the matrix is symmetrical.  Now, I need to pull out the values per shark in matched_data_wide.

#make an empty column for numerical values
genus_matched_data_wide$sor_dist<-NA_real_

for (shark_i in 1:nrow(genus_matched_data_wide)){
  stomach_id<-genus_matched_data_wide$stomach[shark_i]
  fecal_id<-genus_matched_data_wide$fecal[shark_i]
  genus_matched_data_wide$sor_dist[shark_i]<-genus_pairwise_sorensen_distances[stomach_id,fecal_id]
}

# Plot Similarity as a scatter plot
genus_similarity_plot <- ggplot(data = genus_matched_data_wide, aes(x = sample_id, y = sor_dist, fill = sor_dist)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +  # Adjust color scale
  labs(title = "Sørensen Pairwise Similarity between Prey Genuses found in Stomach and Fecal Samples",
       x = "Shark_ID",
       y = "Sorensen Distances") +
  theme_minimal() +
  theme(text = element_text(size = 10),  # Adjust text size
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"))  # Adjust title size and style

print(genus_similarity_plot)


#Do the same Sørensen distance matrix and scatterplot with family tax_glom

#merge motus by prey family
NC_Shark_prey_family_merge.ps <- tax_glom(NC_Sharks.ps, taxrank = "family", NArm = FALSE)
plot_bar(NC_Shark_prey_family_merge.ps, fill = "family") + 
  geom_bar(aes(color=family, fill=family), stat="identity", position="stack")

NC_Sharks.ps
NC_Shark_prey_species_merge.ps
NC_Shark_prey_genus_merge.ps
NC_Shark_prey_family_merge.ps


NC_Shark_prey_family_merge.ps
sample_variables(NC_Shark_prey_family_merge.ps)
rank_names(NC_Shark_prey_family_merge.ps)
#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(NC_Shark_prey_family_merge.ps))
standf = function(x, t=total) round(t * (x / sum(x)))
NC_Sharks.ps.family.med.norm = transform_sample_counts(NC_Shark_prey_family_merge.ps, standf)

NC_Sharks.ps.family.med.norm

#Get just the fecal and stomach paired samples

NC_Shark_prey_family_merge_med.norm.paired.ps<-subset_samples(NC_Sharks.ps.family.med.norm,Paired=="Yes")
NC_Shark_prey_family_merge_med.norm.paired.ps
#re-do Sorensen distances for the family phyloseq object

family_matched_fs_and_sts<- sample_data(NC_Shark_prey_family_merge_med.norm.paired.ps)[,c("sample_id","material")]

family_matched_fs_and_sts$id<- rownames(sample_data(NC_Shark_prey_family_merge_med.norm.paired.ps))
family_matched_fs_and_sts<-data.frame(family_matched_fs_and_sts)

family_matched_data_wide<-family_matched_fs_and_sts%>%
  pivot_wider(id_cols = sample_id,names_from = material, values_from = id)
#View(family_matched_data_wide)


# Extract OTU table data for stomach and fecal samples
family_stomach_data <- otu_table(NC_Shark_prey_family_merge_med.norm.paired.ps)[, family_matched_data_wide$stomach]
family_fecal_data <- otu_table(NC_Shark_prey_family_merge_med.norm.paired.ps)[, family_matched_data_wide$fecal]

dim(family_stomach_data)
dim(family_fecal_data)

colnames(family_stomach_data)
colnames(family_fecal_data)

#calculate distances and put them in a symmetrical matrix
family_pairwise_sorensen_distances<-as.matrix(print(phyloseq::distance(NC_Shark_prey_family_merge_med.norm.paired.ps, method="sor", binary=TRUE, type= "samples", upper=TRUE, diag=TRUE)))

View(family_pairwise_sorensen_distances)
class(family_pairwise_sorensen_distances)

#Yay, the matrix is symmetrical.  Now, I need to pull out the values per shark in matched_data_wide.

#make an empty column for numerical values
family_matched_data_wide$sor_dist<-NA_real_

for (shark_i in 1:nrow(family_matched_data_wide)){
  stomach_id<-family_matched_data_wide$stomach[shark_i]
  fecal_id<-family_matched_data_wide$fecal[shark_i]
  family_matched_data_wide$sor_dist[shark_i]<-family_pairwise_sorensen_distances[stomach_id,fecal_id]
}

# Plot Similarity as a scatter plot
family_similarity_plot <- ggplot(data = family_matched_data_wide, aes(x = sample_id, y = sor_dist, fill = sor_dist)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +  # Adjust color scale
  labs(title = "Sørensen Pairwise Similarity between Prey Families found in Stomach and Fecal Samples",
       x = "Shark_ID",
       y = "Sorensen Distances") +
  theme_minimal() +
  theme(text = element_text(size = 10),  # Adjust text size
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"))  # Adjust title size and style

print(family_similarity_plot)

#Do the same Sørensen distance matrix and scatterplot with order tax_glom

#merge motus by prey order
NC_Shark_prey_order_merge.ps <- tax_glom(NC_Sharks.ps, taxrank = "order", NArm = FALSE)
plot_bar(NC_Shark_prey_order_merge.ps, fill = "order") + 
  geom_bar(aes(color=order, fill=order), stat="identity", position="stack")

NC_Sharks.ps
NC_Shark_prey_species_merge.ps
NC_Shark_prey_genus_merge.ps
NC_Shark_prey_family_merge.ps
NC_Shark_prey_order_merge.ps



NC_Shark_prey_order_merge.ps
sample_variables(NC_Shark_prey_order_merge.ps)
rank_names(NC_Shark_prey_order_merge.ps)
#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(NC_Shark_prey_order_merge.ps))
standf = function(x, t=total) round(t * (x / sum(x)))
NC_Sharks.ps.order.med.norm = transform_sample_counts(NC_Shark_prey_order_merge.ps, standf)

NC_Sharks.ps.order.med.norm

#Get just the fecal and stomach paired samples

#NC_Shark_prey_order_merge_med.norm.paired.ps<-subset_samples(NC_Sharks.ps.order.med.norm,Paired=="Yes")
NC_Shark_prey_order_merge_med.norm.paired.ps<-NC_Sharks.ps.order.med.norm
#re-do Sorensen distances for the order phyloseq object

order_matched_fs_and_sts<- sample_data(NC_Shark_prey_order_merge_med.norm.paired.ps)[,c("sample_id","material")]

order_matched_fs_and_sts$id<- rownames(sample_data(NC_Shark_prey_order_merge_med.norm.paired.ps))
order_matched_fs_and_sts<-data.frame(order_matched_fs_and_sts)

order_matched_data_wide<-order_matched_fs_and_sts%>%
  pivot_wider(id_cols = sample_id,names_from = material, values_from = id)
#View(order_matched_data_wide)

# Extract OTU table data for stomach and fecal samples
order_stomach_data <- otu_table(NC_Shark_prey_order_merge_med.norm.paired.ps)[, order_matched_data_wide$stomach]
order_fecal_data <- otu_table(NC_Shark_prey_order_merge_med.norm.paired.ps)[, order_matched_data_wide$fecal]

dim(order_stomach_data)
dim(order_fecal_data)

colnames(order_stomach_data)
colnames(order_fecal_data)

#calculate distances and put them in a symmetrical matrix
order_pairwise_sorensen_distances<-as.matrix(print(phyloseq::distance(NC_Shark_prey_order_merge_med.norm.paired.ps, method="sor", binary=TRUE, type= "samples", upper=TRUE, diag=TRUE)))

View(order_pairwise_sorensen_distances)
class(order_pairwise_sorensen_distances)

#Yay, the matrix is symmetrical.  Now, I need to pull out the values per shark in matched_data_wide.

#make an empty column for numerical values
order_matched_data_wide$sor_dist<-NA_real_

for (shark_i in 1:nrow(order_matched_data_wide)){
  stomach_id<-order_matched_data_wide$stomach[shark_i]
  fecal_id<-order_matched_data_wide$fecal[shark_i]
  order_matched_data_wide$sor_dist[shark_i]<-order_pairwise_sorensen_distances[stomach_id,fecal_id]
}

# Plot Similarity as a scatter plot
order_similarity_plot <- ggplot(data = order_matched_data_wide, aes(x = sample_id, y = sor_dist, fill = sor_dist)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +  # Adjust color scale
  labs(title = "Sørensen Pairwise Similarity between Prey Orders found in Stomach and Fecal Samples",
       x = "Shark_ID",
       y = "Sorensen Distances") +
  theme_minimal() +
  theme(text = element_text(size = 10),  # Adjust text size
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"))  # Adjust title size and style

print(order_similarity_plot)

#End of Sorensen exercise


# play with ordinations interactively
#install.packages("ggraph") # for taxatree_plots()
#install.packages("DT") # for tax_fix_interactive()
#install.packages("corncob") # for example datasets and beta binomial models
library(microViz)
library(ggraph)
library(DT)
library(corncob)


rank_names(NC_Shark_prey_species_merge_med.norm.paired.ps)
#uncomment the next line for interactive taxonomy fixing to give you new code to paste in
#tax_fix_interactive(NC_Shark_prey_species_merge_med.norm.paired.ps)
NC_Shark_prey_species_merge_med.norm.paired.ps %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )

pseq <- NC_Shark_prey_species_merge_med.norm.paired.ps %>%
  phyloseq_validate()
ord_explore(pseq) 
#code from the above program... 

pseq %>%
  tax_transform(rank = "family", trans = "binary") %>%
  dist_calc(dist = "sor") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "material", fill = "material",
    shape = "circle", alpha = 0.5,
    size = 2
  ) +
  stat_chull(
    ggplot2::aes(colour = material)
  )

####### end of Sorensen exercise ######





#the figure you wanted..
NC_Shark_prey_species_merge_paired.ps.prop <- transform_sample_counts(NC_Shark_prey_species_merge_paired.ps, function(species) species/sum(species))

plot_bar(NC_Shark_prey_species_merge_paired.ps.prop,x="sample_id", fill="genus") + 
  facet_grid(~material~Species..scientific.name., scales="free") +
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")
ggsave("data/Paired/07_Phyloseq/paired_samples_read_proportions_by_shark_species_prey_genus.jpg")


plot_heatmap(NC_Shark_prey_species_merge.ps.prop, method = "NMDS", distance = "bray", taxa.label = "species", sample.label="sample_id")
ggsave("data/Paired/07_Phyloseq/heatmap_sharkID_by_species_proportions.jpg")


plot_bar(NC_Shark_prey_species_merge.ps.prop, x="genus", fill = "genus", facet_grid = Species..scientific.name.~material) +
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

plot_bar(NC_Shark_prey_species_merge.ps.prop, fill="genus") + facet_wrap(~Species..scientific.name., scales= "free_x", nrow=1) +
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

ggsave("data/Paired/07_Phyloseq/read_proportions_by_shark_species_prey_genus.jpg")

plot_heatmap(NC_Shark_prey_species_merge.ps.prop, method = "NMDS", distance = "bray", taxa.label = "genus" )
ggsave("data/Paired/07_Phyloseq/heatmap_sample_by_prey_genus_proportions.jpg")

plot_richness(NC_Shark_prey_species_merge.ps, x="Fork.Length..FL..mm..", color="Species..scientific.name.", measures=c("Observed"))

ggsave("data/Paired/07_Phyloseq/alpha_diversity_prey_richness_by_fork_length_all_species.jpg")

plot_richness(NC_Shark_prey_species_merge.ps, x="Fork.Length..FL..mm..", color="Species..scientific.name.", measures=c("Observed", "Shannon")) + geom_boxplot()

plot_richness(NC_Shark_prey_species_merge.ps, x="Species..scientific.name.", color="Species..scientific.name.", measures=c("Observed", "Shannon")) + geom_boxplot()
ggsave("data/Paired/07_Phyloseq/prey_species_diversity_obs_Shannon_by_Shark_species.jpg")

#following https://alexiscarter.github.io/metab/Phyloseq_script.html

sample_data(NC_Sharks.ps)$Species = factor(sample_data(NC_Sharks.ps)$Species..scientific.name., levels = c('Rhizoprionodon terraenovae', 'Carcharhinus limbatus', 'Carcharhinus acronotus','Sphyrna tiburo'))
sample_data(NC_Sharks.ps)$type = factor(sample_data(NC_Sharks.ps)$material, levels = c('fecal', 'stomach'))

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


NC_Shark_prey_species_merge.ps.1<-tax_glom(NC_Sharks.ps, taxrank = "species", NArm = TRUE)
rank_names(NC_Shark_prey_species_merge.ps.1)
#tax_fix_interactive(NC_Shark_prey_species_merge.ps.1)
NC_Shark_prey_species_merge.ps.1 %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = "", anon_unique = TRUE,
    suffix_rank = "classified"
  )
top20species1 <- names(sort(taxa_sums(NC_Shark_prey_species_merge.ps.1), decreasing=TRUE))[0:20]
ps.top20species1 <- transform_sample_counts(NC_Shark_prey_species_merge.ps.1, function(OTU) OTU/sum(OTU))
ps.top20species1 <- prune_taxa(top20species1, ps.top20species1)
plot_bar(ps.top20species1, x="Species..scientific.name.",fill="species")+ 
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")

ggsave("data/Paired/07_Phyloseq/stacked_bar_chart_top20species_level_prey_by_shark_species.jpg")

NC_Shark_prey_species_merge.ps.1
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 43 taxa and 131 samples ]
#sample_data() Sample Data:       [ 131 samples by 61 sample variables ]
#tax_table()   Taxonomy Table:    [ 43 taxa by 7 taxonomic ranks ]

top43species1 <- names(sort(taxa_sums(NC_Shark_prey_species_merge.ps.1), decreasing=TRUE))[0:43]
ps.top43species1 <- transform_sample_counts(NC_Shark_prey_species_merge.ps.1, function(OTU) OTU/sum(OTU))
ps.top20species1 <- prune_taxa(top43species1, ps.top43species1)
plot_bar(ps.top43species1, x="Species..scientific.name.",fill="species")+ 
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")

#too many species for the legend to make sense

# play with ordinations interactively
pseq <- NC_Shark_prey_species_merge.ps %>%
  phyloseq_validate()
#uncomment the next line for a web browser based way to explore ordinations with these data!
ord_explore(pseq) 
#code from the above program... 

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

ggsave("data/Paired/07_Phyloseq/presence_absence_diet_PCA.jpg")

#compositional transformation of data, showing arrows of the differentiating genuses and the strength of their influence on the PCA plot
pseq %>%
  tax_transform(rank = "genus", trans = "compositional") %>%
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

ggsave("data/Paired/07_Phyloseq/PCA_of_diet_composition_with_arrows.jpg")

#PCoA of Jensen Shannon distance ordination of compositional-transformed data
pseq %>%
  tax_transform(rank = "genus", trans = "compositional") %>%
  dist_calc(dist = "jsd") %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Species..scientific.name.", fill = "Species..scientific.name.",
    shape = "circle", alpha = 0.5,
    size = 2
  ) +
  ggplot2::stat_ellipse(
    ggplot2::aes(colour = Species..scientific.name.)
  )

