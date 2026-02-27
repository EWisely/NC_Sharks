#LULU curation, Metabar filtering, removing the donor species from the results, Combining BerryCrust and MiFish results, and some preliminary Phyloseq visualizations.
#Author: Eldridge Wisely

#install.packages("devtools")
library(devtools)
#install_github("tobiasgf/lulu")
require(lulu)
library(phyloseq)


#setwd("~/Desktop/Savannah/following_obitools3_method/MiFish_pureobi3_lulu_metabar") #deprecated in favor of RProject.

#### Work in an Rproject in a project directory with a directory called data/04_LULU_inputs with the results of script 03_prep_for_LULU.sh inside it.  (i.e. NC-Sharks-diet/data/04_LULU_inputs, NC-Sharks-diet/NC_Sharks.Rproj, NC-Sharks-diet/scripts/04_lulu_metabaR_BerryCrust_NC_SHARKS.R)

load("All_NC_Sharks_Data/data/06_MetabaR_results/After_BerryCrust_obi3_lulu_metabaR.RData")


#First, lulu curation of obi3 results file.

#(not matching merged obi3 sequences back to the results for the lulu step, just results-results for matchlist and results otutable)
otutab <- read.csv("data/04_LULU_inputs/MiFish_named_tab_LULU.txt", sep='\t', header=TRUE, as.is=TRUE, row.names = 1) #was made using the obi3 merged reads matched against the obi3 results 
matchlist <- read.table("data/04_LULU_inputs/MiFish_named_matchlist.txt", header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)
curated_result <- lulu(otutab, matchlist)

curated_result$curated_table # Curated OTU table
curated_result$curated_count # Number of OTUs retained
curated_result$curated_otus # IDs of curated OTUs

curated_result$discarded_count # OTUs discarded
curated_result$otu_map # total - total read count, spread - the number of samples the OTU is present in
# parent_id - ID of OTU with which this OTU was merged (or self)
# curated - ("parent" or "merged"), was this OTU kept as a valid OTU (parent) or merged with another
# rank - The rank of the OTU in terms of decreasing spread and read count

curated_result$original_table # Original OTU table
write.csv(curated_result$curated_table,"data/05_MetabaR_inputs/lulu_curated_MiFish_named.tab")


#Phyloseq section
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
#import data:

#otus
MiFish_lulu_table<-otu_table(curated_result$curated_table, taxa_are_rows = TRUE)
Metabar_formatted_samples<-read.csv("data/05_MetabaR_inputs/MetabaR_formatted_samples_NCSharks.txt", sep="\t")

#samples
Metabar_formatted_pcrs<-read.csv("data/05_MetabaR_inputs/MetabaR_formatted_pcrs_NCSharks_MiFish.txt", sep="\t")
samples.df<-full_join(Metabar_formatted_pcrs,Metabar_formatted_samples,by="sample_id")

MiFish_sample_data<-sample_data(samples.df)

#taxonomy
MiFish_named_obi3_results<-read.csv("data/05_MetabaR_inputs/MiFish_named.tab", sep="\t")
MiFish_named_obi3_results<-MiFish_named_obi3_results%>% select(id, COUNT, BEST_IDENTITY, TAXID, SCIENTIFIC_NAME, ID_STATUS, sequence)

#install.packages("taxonomizr")
library(taxonomizr)
prepareDatabase(getAccessions=FALSE)
#make a taxa table for all obi3 result otus
taxaId<-MiFish_named_obi3_results$TAXID
taxa<-getTaxonomy(taxaId,'nameNode.sqlite')
print(taxa)
class(taxa)
taxa.df<-as.data.frame.array(taxa)
rownames(taxa.df)<-MiFish_named_obi3_results$id

#make a taxa table for only the lulu-curated otus
MiFish_lulu_table.df<-as.data.frame(MiFish_lulu_table)
MiFish_lulu_table.df<-rownames_to_column(MiFish_lulu_table.df, var ="id")
taxa.df<-rownames_to_column(taxa.df, var="id")
taxa.df.lulu<- semi_join(taxa.df,MiFish_lulu_table.df,by="id")

#see which taxa were discarded by lulu
#taxa.df.discards<-anti_join(taxa.df,MiFish_lulu_table.df,by="id")

#do something about the negatives

#Metabar!

#format the motus
#taxa.df.lulu has the right motus and taxonomic levels, MiFish_named_obi3_motus has the rest of the info I want but too many motus join by ="id"

Metabar_formatted_motus<- inner_join(MiFish_named_obi3_results,taxa.df.lulu,by="id")
Metabar_formatted_motus<-column_to_rownames(Metabar_formatted_motus, var="id")


#flip rows and columns (transpose) of MiFish_lulu_table.df to make the reads table
MiFish_lulu_table.df<-as.data.frame(MiFish_lulu_table)

Metabar_formatted_reads<-t(MiFish_lulu_table.df)

Metabar_formatted_reads<-as.data.frame(Metabar_formatted_reads)

Metabar_formatted_reads<-as.matrix(Metabar_formatted_reads)

Metabar_formatted_samples<-column_to_rownames(Metabar_formatted_samples, var="sample_id")

Metabar_formatted_pcrs<-column_to_rownames(Metabar_formatted_pcrs, var="sample")



#NC_MiFish<-metabarlist_generator(reads=Metabar_formatted_reads, Metabar_formatted_motus, pcrs=Metabar_formatted_pcrs, samples=Metabar_formatted_samples)
#this seems to work, but doesn't identify any pcr contaminants which makes me distrust it, so, I'm going to export to csv and look at them in excel and reimport with tabfilestometabarlist


write.csv(Metabar_formatted_reads, "test-MiFish_Metabar_formatted_reads_pureobi3_lulu.csv")
write.csv(Metabar_formatted_motus, "test-MiFish_Metabar_formatted_motus_pureobi3_lulu.csv")
write.csv(Metabar_formatted_pcrs, "test-MiFish_Metabar_formatted_pcrs_pureobi3_lulu.csv")
write.csv(Metabar_formatted_samples, "test-MiFish_Metabar_formatted_samples_pureobi3_lulu.csv")


#samples file has lots of extra periods in it and the first line doesn't say sample_id, so I saved MetabaR_formatted_samples_NCSharks.txt to this folder to use instead.
#pcrs file looked exactly the same as a previous one I used for MiFish analysis except the first column label (sample) was gone so I put it back and saved it as a tab delimited text. MiFish_Metabar_formatted_pcrs_pureobi3_lulu.txt
#pcrs had St9 and reads did not.  So I removed the St9 row from pcrs and saved as MiFish_Metabar_formatted_pcrs_pureobi3_lulu_noSt9.txt
#motus needed some rearranging to make it like one of the MiFish ones that worked in the past.  #I put the count field second and made it "count" instead of "COUNT", then the next field is seq_length and I calculated that in excel from the sequence field, which I put at the end of the taxonomy fields. Saved as MiFish_Metabar_formatted_motus_pureobi3_lulu.txt
#reads says St#.St# and pcrs says St#-St#!  changed reads to "-"instead of ".". saved as MiFish_Metabar_formatted_reads_pureobi3_lulu.txt


NC_MiFish<-tabfiles_to_metabarlist(
  file_reads= "data/05_MetabaR_inputs/MiFish_Metabar_formatted_reads_pureobi3_lulu.txt",
  file_motus= "data/05_MetabaR_inputs/MiFish_Metabar_formatted_motus_pureobi3_lulu.txt",
  file_pcrs ="data/05_MetabaR_inputs/MiFish_Metabar_formatted_pcrs_pureobi3_lulu_noSt9.txt",
  file_samples ="data/05_MetabaR_inputs/MetabaR_formatted_samples_NCSharks.txt",
  files_sep = "\t"
)

#same error so I copied the pcrs entries and pasted into a new file saved under the same name in case there was a weird empty row or something.
#still same error
#copied and pasted reads even though I hadn't messed with that.
#same error
#deleted "sample" from the first rowname of each, and still got the same error.
#aha! reads says St#.St# and pcrs says St#-St#!  changed reads to "-"instead of ".".

#yay! got rid of that error!

#now, Error in check_metabarlist(out) : 
#table `motus` in out has empty column names

#I added id to the first column name of the otus.
#still the same problem.
#deleted the empty COUNT column that I had cut and pasted the contents of.

#fixed! Finally no errors!



#Arguments
#reads	
#MOTU abundance table. Rows and rownames of the table should correspond to PCRs and their names respectively. Columns and colnames should correspond to MOTUs and their names. Rownames in this table should correspond to PCR names respectively.

#motus	
#MOTU characteristics table (e.g. taxonomy, sequence, etc.). Rows and rownames of the table should correspond to MOTUs and their names respectively, and the columns to their characteristics. Mandatory fields: 'sequence', i.e. the sequence representative of the MOTU.

#pcrs	
#PCR characteristics table (e.g. tags, primers, plate wells, etc.). Rows and rownames of the table should correspond to PCRs and their names respectively, and the columns to their characteristics. Mandatory fields: (i) 'sample_id', i.e. the name of each biological sample. (ii) 'type', i.e. the type of PCR; can be 'sample' or 'control'. (iii) 'control_type', i.e. the type of control if applicable. Should be either: 'NA' for samples, 'extraction' for extraction negative controls, 'pcr' for PCR negative controls, 'sequencing' for sequencing negative controls (e.g. unused tag combinations), or 'positive' for positive controls.

#samples	
#Sample characteristics table. Rows and rownames of the table should correspond to biological samples and their names respectively, and the columns to their environnemental characteristics.


# Compute the number of reads per pcr
NC_MiFish$pcrs$nb_reads <- rowSums(NC_MiFish$reads)

# Compute the number of motus per pcr
NC_MiFish$pcrs$nb_motus <- rowSums(NC_MiFish$reads>0)

# Load requested package for plotting
library(ggplot2)
library(reshape2)
library(tidyverse)

# Create an input table (named check1) for ggplot of 3 columns: 
#  (i) control type 
#  (ii) a vector indicated whether it corresponds to nb_reads or nb_motus, 
#  (iii) the corresponding values.

check1 <- reshape2::melt(NC_MiFish$pcrs[,c("control_type", "nb_reads", "nb_motus")])

ggplot(data <- check1, aes(x=control_type, y=value, color=control_type)) + 
  geom_boxplot() + theme_bw() + 
  geom_jitter(alpha=0.2) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable, scales = "free_y") + 
  theme(axis.text.x = element_text(angle=45, h=1))

# Using the nb_reads and nb_motus defined previously in the NC_MiFish$pcrs table

ggplot(NC_MiFish$pcrs, aes(x=nb_reads, y=nb_motus, color = control_type)) + 
  geom_point() + theme_bw() + 
  scale_y_log10() + scale_x_log10() + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey")

ggpcrplate(NC_MiFish, FUN = function(m){rowSums(m$reads)}, legend_title = "# of reads per PCR")


# Here the list of all tag/indices used in the experiment 
# is available in the column "tag_rev" of the NC_MiFish$pcrs table
tag.list <- as.character(unique(NC_MiFish$pcrs$tag_fwd))
ggpcrtag(NC_MiFish, 
         legend_title = "# of reads per PCR", 
         FUN = function(m) {rowSums(m$reads)},
         taglist = tag.list) 


NC_MiFish.raref = hill_rarefaction(NC_MiFish, nboot = 20, nsteps = 10)
head(NC_MiFish.raref$hill_table)
gghill_rarefaction(NC_MiFish.raref) 


# Define a vector containing the Material info for each pcrs 
material <- NC_MiFish$pcrs$material

# Use of gghill_rarefaction requires a vector with named pcrs
material <- setNames(material,rownames(NC_MiFish$pcrs))

# Plot
p <- gghill_rarefaction(NC_MiFish.raref, group=material)
p + scale_fill_manual(values = c("goldenrod4", "brown4", "grey")) +
  scale_color_manual(values = c("goldenrod4", "brown4", "grey")) +
  labs(color="Material type")


# Identifying extraction contaminants(changed to pcr since that's the only control type we have)
NC_MiFish <- contaslayer(NC_MiFish, 
                             control_types = "pcr",
                             output_col = "not_a_pcr_conta")

table(NC_MiFish$motus$not_a_pcr_conta)
#I blasted the top 3 contaminants and they came back with very low confidence in catsharks and zebrasharks.

# Identify the most common contaminant
# get contaminant ids
id <- !NC_MiFish$motus$not_a_pcr_conta
max.conta <- rownames(NC_MiFish$motus[id,])[which.max(NC_MiFish$motus[id, "count"])]

#... and its distribution and relative abundance in each pcr
ggpcrplate(NC_MiFish, legend_title = "#reads of most \nabundant contaminant",
           FUN = function(m) {m$reads[, max.conta]/rowSums(m$reads)})

# Compute relative abundance of all pcr contaminants together 
a <- data.frame(conta.relab = rowSums(NC_MiFish$reads[,!NC_MiFish$motus$not_a_pcr_conta]) / 
                  rowSums(NC_MiFish$reads))
# Add information on control types
a$control_type <- NC_MiFish$pcrs$control_type[match(rownames(a), rownames(NC_MiFish$pcrs))]

ggplot(a, aes(x=control_type, y=conta.relab, color=control_type)) + 
  geom_boxplot() + geom_jitter(alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x=NULL, y="Prop. Reads (log10)") + 
  theme_bw() + 
  scale_y_log10()


# flag pcrs with total contaminant relative abundance > 20% of reads)
NC_MiFish$pcrs$low_contamination_level <- 
  ifelse(a$conta.relab[match(rownames(NC_MiFish$pcrs), rownames(a))]>2e-1,  F, T)

# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(NC_MiFish$pcrs$low_contamination_level) / nrow(NC_MiFish$pcrs)

#next the tutorial excludes non-eukaryotes and plots the similarity scores and removes the ones with low scores.  
#Whatever metabarcoding pipeline hasn't done that yet isn't worth it's salt IMO

#Flag MOTUs corresponding to target (TRUE) vs. non-target (FALSE) taxa 
NC_MiFish$motus$target_taxon <- grepl("TRUE", NC_MiFish$motus$ID_STATUS)#I changed this from Eukaryota and path because I am using the obitools3 results instead of the sintax vsearch results

# Proportion of each of these over total number of MOTUs
table(NC_MiFish$motus$target_taxon) / nrow(NC_MiFish$motus)

# Intersection with extraction contaminant flags (not contaminant = T)
table(NC_MiFish$motus$target_taxon, 
      NC_MiFish$motus$not_a_pcr_conta)


#remove non-target taxa and MF control contaminants
NC_MiFish <- contaslayer(NC_MiFish, 
                             method= "all",
                             control_types = "pcr",
                             controls = c("NCFSMFneg","NCSTMFneg"),
                             output_col = "not_a_MF_pcr_conta")

table(NC_MiFish$motus$not_a_pcr_conta)
#I'd like to be able to look at which 65 were FALSE
table(NC_MiFish$motus$not_a_MF_pcr_conta)
#22 were false
table(NC_MiFish$motus$target_taxon, 
      NC_MiFish$motus$not_a_MF_pcr_conta)
table(NC_MiFish$motus$not_a_MF_pcr_conta, 
      NC_MiFish$motus$not_a_pcr_conta)

#all of the pcr contas are MF_pcr_contas, and none are target taxa! yay!

#remove sequences from our donor shark species

NC_MiFish$motus$Sphyrna_tiburo_motus <- grepl("Sphyrna", NC_MiFish$motus$SCIENTIFIC_NAME)
NC_MiFish$motus$C_acronotus_motus <- grepl("Carcharhinus", NC_MiFish$motus$SCIENTIFIC_NAME)
NC_MiFish$motus$C_limbatus_motus <- grepl("Carcharhinus", NC_MiFish$motus$SCIENTIFIC_NAME)
NC_MiFish$motus$R_terranovae_motus <- grepl("Rhizoprionodon", NC_MiFish$motus$SCIENTIFIC_NAME)

NC_MiFish_no_donors <-subset_metabarlist(NC_MiFish, table ="motus", indices =rowSums(NC_MiFish$motus[,c("Sphyrna_tiburo_motus", "C_acronotus_motus","C_limbatus_motus", "R_terranovae_motus")]) == 0 )

summary_metabarlist(NC_MiFish_no_donors)

NC_MiFish_no_donors1 <- subset_metabarlist(NC_MiFish_no_donors, table = "pcrs",
                                          indices = NC_MiFish_no_donors$pcrs$nb_reads>0)

summary_metabarlist(NC_MiFish_no_donors1)

#332 motus remaining out of 384 from lulu, 158 pcrs (all of them)

summary_metabarlist(NC_MiFish)

NC_MiFish1 <- subset_metabarlist(NC_MiFish_no_donors, "motus", 
                                     indices = rowSums(NC_MiFish_no_donors$motus[,c("not_a_MF_pcr_conta", "target_taxon")]) == 2)


summary_metabarlist(NC_MiFish1)

#138 motus remaining average 11.8

# Plot the unweighted distribution of MOTUs similarity scores 
a <- 
  ggplot(NC_MiFish1$motus, aes(x=BEST_IDENTITY)) + 
  geom_histogram(color="grey", fill="white", bins=20) + 
  geom_vline(xintercept = 0.97, col="orange", lty=2) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(x="% similarity against best match", y="# MOTUs")

# Same for the weighted distribution
b <- 
  ggplot(NC_MiFish1$motus, 
         aes(x=BEST_IDENTITY, y = after_stat(count), weight = count)) + 
  geom_histogram(color="grey", fill="white", bins=20) + 
  geom_vline(xintercept = 0.97, col="orange", lty=2) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(x="% similarity against best match", y="# Reads")

# Combine plots into one
library(cowplot)
ggdraw() + 
  draw_plot(a, x=0, y=0, width = 0.5) + 
  draw_plot(b, x=0.5, y=0, width = 0.5)
#there are some motus with less than 99% identity but hardly any reads under 99%.  These could be rare species or PCR errors that lulu didn't catch

########This is where I can play with different values for ID#####

# Flag not degraded (TRUE) vs. potentially degraded sequences (FALSE)
NC_MiFish1$motus$not_degraded <-
  ifelse(NC_MiFish1$motus$BEST_IDENTITY < 0.98, F, T)

# Proportion of each of these over total number of MOTUs
table(NC_MiFish1$motus$not_degraded) / nrow(NC_MiFish1$motus)

# Flag not degraded (TRUE) vs. potentially degraded sequences (FALSE)
NC_MiFish$motus$not_degraded <-
  ifelse(NC_MiFish$motus$BEST_IDENTITY < 0.98, F, T)

# Proportion of each of these over total number of MOTUs
table(NC_MiFish$motus$not_degraded) / nrow(NC_MiFish$motus)
# for NC_MiFish
#> 
#    FALSE      TRUE 
#0.8042813 0.1957187 
#> 
#> For NC_MiFish1
#
#FALSE      TRUE 
#0.1351351 0.8648649 

# Intersection with other flags
table(NC_MiFish$motus$target_taxon, 
      NC_MiFish$motus$not_a_pcr_conta, 
      NC_MiFish$motus$not_degraded)


#       FALSE TRUE
#FALSE    65  129
#TRUE      0    5


#       FALSE TRUE
#FALSE     0    0
#TRUE      0  185

#Moving on to detecting PCR outliers by sequencing depth
ggplot(NC_MiFish$pcrs, aes(nb_reads)) +
  geom_histogram(bins=40, color="grey", fill="white") + 
  geom_vline(xintercept = 2e3, lty=2, color="orange") + # threshold
  scale_x_log10() + 
  labs(x="# Reads (with all MOTUs and PCRs)", 
       y="# PCRs") +
  theme_bw() + 
  theme(panel.grid = element_blank())


#play with the xintercept above to visualize your cutoff level

# Flag pcrs with an acceptable sequencing depth (TRUE) or inacceptable one (FALSE)
NC_MiFish$pcrs$seqdepth_ok <- ifelse(NC_MiFish$pcrs$nb_reads < 2e3, F, T)

# Proportion of each of these over total number of pcrs, control excluded
table(NC_MiFish$pcrs$seqdepth_ok[NC_MiFish$pcrs$type=="sample"]) /
  nrow(NC_MiFish$pcrs[NC_MiFish$pcrs$type=="sample",])
#and for the semi-cleaned dataset
NC_MiFish1$pcrs$seqdepth_ok <- ifelse(NC_MiFish1$pcrs$nb_reads < 2e3, F, T)

# Proportion of each of these over total number of pcrs, control excluded
table(NC_MiFish1$pcrs$seqdepth_ok[NC_MiFish1$pcrs$type=="sample"]) /
  nrow(NC_MiFish1$pcrs[NC_MiFish1$pcrs$type=="sample",])


#They use this next function to compare the similarity of biological controls with the expectation they will be more similar to each other than to other samples.  
#We can use this to check the similarity of fecals to their associated stomachs

# Subsetting the metabarlist
NC_MiFish_sub <- subset_metabarlist(NC_MiFish, 
                                        table = "pcrs", 
                                        indices = NC_MiFish$pcrs$nb_reads>0 & (
                                          is.na(NC_MiFish$pcrs$control_type) |
                                            NC_MiFish$pcrs$control_type=="positive"))

# First visualization
comp1 = pcr_within_between(NC_MiFish_sub)
check_pcr_thresh(comp1)

#cool, this looks good.  Now... for flagging not-well-replicated pcrs.  Probably shouldn't flag these for removal since that's not appropriate for this experimental design

# Subsetting the metabarlist
NC_MiFish1_sub <- subset_metabarlist(NC_MiFish1, 
                                         table = "pcrs", 
                                         indices = NC_MiFish1$pcrs$nb_reads>0 & (
                                           is.na(NC_MiFish$pcrs$control_type) |
                                             NC_MiFish$pcrs$control_type=="positive"))

# First visualization
comp2 = pcr_within_between(NC_MiFish1_sub)
check_pcr_thresh(comp2)

#after removing the non-target , donor species, and MiFish contamination motus, the distance between samples and within samples stayed the same but the density of markers was lower and more even.


#now for the "lowering tag-jumps" section:
# Define a vector of thresholds to test
thresholds <- c(0,1e-4,1e-3, 1e-2, 2e-2,3e-2, 5e-2) 

# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(NC_MiFish,x, method = "substract"))
names(tests) <- paste("t_", thresholds, sep="")

# Format the data for ggplot with amount of reads at each threshold
tmp <- reshape2::melt(as.matrix(do.call("rbind", lapply(tests, function(x) rowSums(x$reads)))))
colnames(tmp) <- c("threshold", "sample", "abundance")

# Add richness in MOTUs at each threshold
tmp$richness <-
  reshape2::melt(as.matrix(do.call("rbind", lapply(tests, function(x) {
    rowSums(x$reads > 0)
  }))))$value

# Add control type information on pcrs and make data curation threshold numeric
tmp$controls <- NC_MiFish$pcrs$control_type[match(tmp$sample, rownames(NC_MiFish$pcrs))]
tmp$threshold <- as.numeric(gsub("t_", "", tmp$threshold))

# New table formatting for ggplot
tmp2 <- reshape2::melt(tmp, id.vars=colnames(tmp)[-grep("abundance|richness", colnames(tmp))])

ggplot(tmp2, aes(x=as.factor(threshold), y=value)) + 
  geom_boxplot(color="grey40") + 
  geom_vline(xintercept = which(levels(as.factor(tmp2$threshold)) == "0.01"), col="orange", lty=2) + 
  geom_jitter(aes(color=controls), width = 0.2, alpha=0.5) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable+controls, scale="free_y", ncol=5) + 
  theme_bw() + 
  scale_y_log10() +
  labs(x="MOTU pcr : total abundance filtering threshold", y="# Reads/MOTUs") + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle=40, h=1), 
        legend.position = "none")
#Removed 16 rows containing non-finite values (`stat_boxplot()`).


#again using the tutorial code which uses cut instead of "substract"

# Define a vector of thresholds to test
thresholds <- c(0,1e-4,1e-3, 1e-2, 2e-2,3e-2, 5e-2) 

# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(NC_MiFish,x))
names(tests) <- paste("t_", thresholds, sep="")

# Format the data for ggplot with amount of reads at each threshold
tmp <- reshape2::melt(as.matrix(do.call("rbind", lapply(tests, function(x) rowSums(x$reads)))))
colnames(tmp) <- c("threshold", "sample", "abundance")

# Add richness in MOTUs at each threshold
tmp$richness <-
  reshape2::melt(as.matrix(do.call("rbind", lapply(tests, function(x) {
    rowSums(x$reads > 0)
  }))))$value

# Add control type information on pcrs and make data curation threshold numeric
tmp$controls <- NC_MiFish$pcrs$control_type[match(tmp$sample, rownames(NC_MiFish$pcrs))]
tmp$threshold <- as.numeric(gsub("t_", "", tmp$threshold))

# New table formatting for ggplot
tmp2 <- reshape2::melt(tmp, id.vars=colnames(tmp)[-grep("abundance|richness", colnames(tmp))])

ggplot(tmp2, aes(x=as.factor(threshold), y=value)) + 
  geom_boxplot(color="grey40") + 
  geom_vline(xintercept = which(levels(as.factor(tmp2$threshold)) == "0.01"), col="orange", lty=2) + 
  geom_jitter(aes(color=controls), width = 0.2, alpha=0.5) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable+controls, scale="free_y", ncol=5) + 
  theme_bw() + 
  scale_y_log10() +
  labs(x="MOTU pcr : total abundance filtering threshold", y="# Reads/MOTUs") + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle=40, h=1), 
        legend.position = "none")

#Removed 10 rows containing non-finite values (`stat_boxplot()`).

#Did this to the semi-cleaned NC_MiFish1 metabarlist to get the BCneglevel of contamination out of the MFneg samples.

# Define a vector of thresholds to test
thresholds <- c(0,1e-4,1e-3, 1e-2, 2e-2,3e-2, 4e-2, 5e-2) 

# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(NC_MiFish1,x))
names(tests) <- paste("t_", thresholds, sep="")

# Format the data for ggplot with amount of reads at each threshold
tmp <- reshape2::melt(as.matrix(do.call("rbind", lapply(tests, function(x) rowSums(x$reads)))))
colnames(tmp) <- c("threshold", "sample", "abundance")

# Add richness in MOTUs at each threshold
tmp$richness <-
  reshape2::melt(as.matrix(do.call("rbind", lapply(tests, function(x) {
    rowSums(x$reads > 0)
  }))))$value

# Add control type information on pcrs and make data curation threshold numeric
tmp$controls <- NC_MiFish1$pcrs$control_type[match(tmp$sample, rownames(NC_MiFish1$pcrs))]
tmp$threshold <- as.numeric(gsub("t_", "", tmp$threshold))

# New table formatting for ggplot
tmp2 <- reshape2::melt(tmp, id.vars=colnames(tmp)[-grep("abundance|richness", colnames(tmp))])

ggplot(tmp2, aes(x=as.factor(threshold), y=value)) + 
  geom_boxplot(color="grey40") + 
  geom_vline(xintercept = which(levels(as.factor(tmp2$threshold)) == "0.01"), col="orange", lty=2) + 
  geom_jitter(aes(color=controls), width = 0.2, alpha=0.5) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable+controls, scale="free_y", ncol=5) + 
  theme_bw() + 
  scale_y_log10() +
  labs(x="MOTU pcr : total abundance filtering threshold", y="# Reads/MOTUs") + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle=40, h=1), 
        legend.position = "none")

#Removed 118 rows containing non-finite values (`stat_boxplot()`).
#after adding the code to remove donor sharks from MiFish1, 
#Removed 638 rows containing non-finite values (`stat_boxplot()`). But still needed the 0.04 threshold for removing all BCneg contaminants.


#Yay, this looks awesome threshold 0.04 removes all BCneg motus and most of the MFneg motus in the cleaned dataset.


#Summarizing noise in the semi-cleaned dataset

# Create a table of MOTUs quality criteria 
# noise is identified as FALSE in soil_euk, the "!" transforms it to TRUE
motus.qual <- !NC_MiFish1$motus[,c("not_a_pcr_conta", "target_taxon", "not_degraded")]
colnames(motus.qual) <- c("pcr_conta", "untargeted_taxon", "degraded_seq")

# Proportion of MOTUs potentially artifactual (TRUE) based on the criteria used
prop.table(table(apply(motus.qual, 1, sum) > 0))

# Corresponding proportion of artifactual reads (TRUE)
prop.table(xtabs(NC_MiFish1$motus$count~apply(motus.qual, 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(motus.qual, 2, sum) / nrow(motus.qual)
apply(motus.qual, 2, function(x) sum(NC_MiFish1$motus$count[x])/sum(NC_MiFish1$motus$count))

tmp.motus <- 
  apply(sapply(1:ncol(motus.qual), function(x) {
    ifelse(motus.qual[,x]==T, colnames(motus.qual)[x], NA)}), 1, function(x) {
      paste(sort(unique(x)), collapse = "|")
    })
tmp.motus <- as.data.frame(gsub("^$", "not_artefactual", tmp.motus))
colnames(tmp.motus) <-  "artefact_type"

ggplot(tmp.motus, aes(x=1, fill=artefact_type)) +
  geom_bar() +  xlim(0, 2) +
  labs(fill="Artifact type") + 
  coord_polar(theta="y") + theme_void() + 
  scale_fill_brewer(palette = "Set3") + 
  theme(legend.direction = "vertical") + 
  ggtitle("MOTUs artefacts overview")


#and the full dataset

#again using the tutorial code which uses cut instead of "substract"

# Define a vector of thresholds to test
thresholds <- c(0,1e-4,1e-3, 1e-2, 2e-2,3e-2, 4e-2,5e-2) 

# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(NC_MiFish,x))
names(tests) <- paste("t_", thresholds, sep="")

# Format the data for ggplot with amount of reads at each threshold
tmp <- reshape2::melt(as.matrix(do.call("rbind", lapply(tests, function(x) rowSums(x$reads)))))
colnames(tmp) <- c("threshold", "sample", "abundance")

# Add richness in MOTUs at each threshold
tmp$richness <-
  reshape2::melt(as.matrix(do.call("rbind", lapply(tests, function(x) {
    rowSums(x$reads > 0)
  }))))$value

# Add control type information on pcrs and make data curation threshold numeric
tmp$controls <- NC_MiFish$pcrs$control_type[match(tmp$sample, rownames(NC_MiFish$pcrs))]
tmp$threshold <- as.numeric(gsub("t_", "", tmp$threshold))

# New table formatting for ggplot
tmp2 <- reshape2::melt(tmp, id.vars=colnames(tmp)[-grep("abundance|richness", colnames(tmp))])

ggplot(tmp2, aes(x=as.factor(threshold), y=value)) + 
  geom_boxplot(color="grey40") + 
  geom_vline(xintercept = which(levels(as.factor(tmp2$threshold)) == "0.01"), col="orange", lty=2) + 
  geom_jitter(aes(color=controls), width = 0.2, alpha=0.5) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable+controls, scale="free_y", ncol=5) + 
  theme_bw() + 
  scale_y_log10() +
  labs(x="MOTU pcr : total abundance filtering threshold", y="# Reads/MOTUs") + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle=40, h=1), 
        legend.position = "none")

#Removed 14 rows containing non-finite values (`stat_boxplot()`).


# Create a table of MOTUs quality criteria for the semi-cleaned dataset
# noise is identified as FALSE in NC_MiFish, the "!" transforms it to TRUE
motus.qual <- !NC_MiFish1$motus[,c("not_a_pcr_conta", "target_taxon", "not_degraded")]
colnames(motus.qual) <- c("pcr_conta", "untargeted_taxon", "degraded_seq")

# Proportion of MOTUs potentially artifactual (TRUE) based on the criteria used
prop.table(table(apply(motus.qual, 1, sum) > 0))

# Corresponding proportion of artifactual reads (TRUE)
prop.table(xtabs(NC_MiFish$motus$count~apply(motus.qual, 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(motus.qual, 2, sum) / nrow(motus.qual)
apply(motus.qual, 2, function(x) sum(NC_MiFish$motus$count[x])/sum(NC_MiFish$motus$count))

tmp.motus <- 
  apply(sapply(1:ncol(motus.qual), function(x) {
    ifelse(motus.qual[,x]==T, colnames(motus.qual)[x], NA)}), 1, function(x) {
      paste(sort(unique(x)), collapse = "|")
    })
tmp.motus <- as.data.frame(gsub("^$", "not_artefactual", tmp.motus))
colnames(tmp.motus) <-  "artefact_type"

ggplot(tmp.motus, aes(x=1, fill=artefact_type)) +
  geom_bar() +  xlim(0, 2) +
  labs(fill="Artifact type") + 
  coord_polar(theta="y") + theme_void() + 
  scale_fill_brewer(palette = "Set3") + 
  theme(legend.direction = "vertical") + 
  ggtitle("MOTUs artefacts overview")

#SAME FOR PCRS

#after cleaning out the non-target taxa and MF neg conta

# Create a table of pcrs quality criteria 
# noise is identified as FALSE in NC_MiFish1, the "!" transforms it to TRUE
pcrs.qual <- !NC_MiFish1$pcrs[,c("low_contamination_level", "seqdepth_ok")]#, "replicating_pcr")]
colnames(pcrs.qual) <- c("high_contamination_level", "low_seqdepth")#, "outliers")

# Proportion of pcrs potentially artifactual (TRUE) based on the criteria used
# excluding controls
prop.table(table(apply(pcrs.qual[NC_MiFish1$pcrs$type=="sample",], 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(pcrs.qual[NC_MiFish1$pcrs$type=="sample",], 2, sum) / nrow(pcrs.qual[NC_MiFish1$pcrs$type=="sample",])

tmp.pcrs <- 
  apply(sapply(1:ncol(pcrs.qual), function(x) {
    ifelse(pcrs.qual[NC_MiFish1$pcrs$type=="sample",x]==T, 
           colnames(pcrs.qual)[x], NA)}), 1, function(x) {
             paste(sort(unique(x)), collapse = "|")
           })
tmp.pcrs <- as.data.frame(gsub("^$", "not_artefactual", tmp.pcrs))

colnames(tmp.pcrs) <- "artefact_type"

ggplot(tmp.pcrs, aes(x=1, fill=artefact_type)) +
  geom_bar() +  xlim(0, 2) +
  labs(fill="Artifact type") + 
  coord_polar(theta="y") + theme_void() + 
  scale_fill_brewer(palette = "Set3") + 
  theme(legend.direction = "vertical") +
  ggtitle("PCR artefacts overview")


#for the full dataset

# Create a table of MOTUs quality criteria 
# noise is identified as FALSE in NC_MiFish, the "!" transforms it to TRUE
motus.qual <- !NC_MiFish$motus[,c("not_a_pcr_conta", "target_taxon", "not_degraded")]
colnames(motus.qual) <- c("pcr_conta", "untargeted_taxon", "degraded_seq")

# Proportion of MOTUs potentially artifactual (TRUE) based on the criteria used
prop.table(table(apply(motus.qual, 1, sum) > 0))

# Corresponding proportion of artifactual reads (TRUE)
prop.table(xtabs(NC_MiFish$motus$count~apply(motus.qual, 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(motus.qual, 2, sum) / nrow(motus.qual)
apply(motus.qual, 2, function(x) sum(NC_MiFish$motus$count[x])/sum(NC_MiFish$motus$count))

tmp.motus <- 
  apply(sapply(1:ncol(motus.qual), function(x) {
    ifelse(motus.qual[,x]==T, colnames(motus.qual)[x], NA)}), 1, function(x) {
      paste(sort(unique(x)), collapse = "|")
    })
tmp.motus <- as.data.frame(gsub("^$", "not_artefactual", tmp.motus))
colnames(tmp.motus) <-  "artefact_type"

ggplot(tmp.motus, aes(x=1, fill=artefact_type)) +
  geom_bar() +  xlim(0, 2) +
  labs(fill="Artifact type") + 
  coord_polar(theta="y") + theme_void() + 
  scale_fill_brewer(palette = "Set3") + 
  theme(legend.direction = "vertical") + 
  ggtitle("MOTUs artefacts overview")

# Create a table of pcrs quality criteria 
# noise is identified as FALSE in NC_MiFish, the "!" transforms it to TRUE
pcrs.qual <- !NC_MiFish$pcrs[,c("low_contamination_level", "seqdepth_ok")]#, "replicating_pcr")]
colnames(pcrs.qual) <- c("high_contamination_level", "low_seqdepth")#, "outliers")

# Proportion of pcrs potentially artifactual (TRUE) based on the criteria used
# excluding controls
prop.table(table(apply(pcrs.qual[NC_MiFish$pcrs$type=="sample",], 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(pcrs.qual[NC_MiFish$pcrs$type=="sample",], 2, sum) / nrow(pcrs.qual[NC_MiFish$pcrs$type=="sample",])

tmp.pcrs <- 
  apply(sapply(1:ncol(pcrs.qual), function(x) {
    ifelse(pcrs.qual[NC_MiFish$pcrs$type=="sample",x]==T, 
           colnames(pcrs.qual)[x], NA)}), 1, function(x) {
             paste(sort(unique(x)), collapse = "|")
           })
tmp.pcrs <- as.data.frame(gsub("^$", "not_artefactual", tmp.pcrs))

colnames(tmp.pcrs) <- "artefact_type"

ggplot(tmp.pcrs, aes(x=1, fill=artefact_type)) +
  geom_bar() +  xlim(0, 2) +
  labs(fill="Artifact type") + 
  coord_polar(theta="y") + theme_void() + 
  scale_fill_brewer(palette = "Set3") + 
  theme(legend.direction = "vertical") +
  ggtitle("PCR artefacts overview")




#remove artefacts and controls

# Use tag-jump corrected metabarlist with the threshold identified above
tmp <- tests[["t_0.04"]]

# Subset on MOTUs: we keep motus that are defined as TRUE following the 
# three criteria below (sum of three TRUE is equal to 3 with the rowSums function)
tmp <- subset_metabarlist(tmp, "motus", 
                          indices = rowSums(tmp$motus[,c("not_a_pcr_conta", "target_taxon",
"not_degraded")]) == 3)
summary_metabarlist(tmp)

#185 motus remaining not using the 98% similarity cutoff.

# Subset on pcrs and keep no controls 
NC_MiFish_clean <- subset_metabarlist(tmp, "pcrs", 
                                          indices = rowSums(tmp$pcrs[,c("low_contamination_level", "seqdepth_ok")]) == 2 & 
                                            tmp$pcrs$type == "sample")
summary_metabarlist(NC_MiFish_clean)

#182 motus remaining

if(sum(colSums(NC_MiFish_clean$reads)==0)>0){print("empty motus present")}
if(sum(colSums(NC_MiFish_clean$reads)==0)>0){print("empty pcrs present")}

NC_MiFish_clean$motus$count = colSums(NC_MiFish_clean$reads)
NC_MiFish_clean$pcrs$nb_reads_postmetabaR = rowSums(NC_MiFish_clean$reads)
NC_MiFish_clean$pcrs$nb_motus_postmetabaR = rowSums(ifelse(NC_MiFish_clean$reads>0, T, F))

check <- reshape2::melt(NC_MiFish_clean$pcrs[,c("nb_reads", "nb_reads_postmetabaR", 
                                                    "nb_motus", "nb_motus_postmetabaR")])
check$type <- ifelse(grepl("motus", check$variable), "richness", "abundance")

ggplot(data = check, aes(x = variable, y = value)) +
  geom_boxplot( color = "darkgrey") +
  geom_jitter(alpha=0.1, color = "darkgrey") +
  theme_bw() +
  facet_wrap(~type, scales = "free", ncol = 5) +
  theme(axis.text.x = element_text(angle=45, h=1))


# Compute the number of reads per pcr
NC_MiFish_clean$pcrs$nb_reads <- rowSums(NC_MiFish_clean$reads)

# Compute the number of motus per pcr
NC_MiFish_clean$pcrs$nb_motus <- rowSums(NC_MiFish_clean$reads>0)

# Load requested package for plotting
library(ggplot2)
library(reshape2)
library(tidyverse)

# Create an input table (named check2) for ggplot of 3 columns: 
#  (i) control type 
#  (ii) a vector indicated whether it corresponds to nb_reads or nb_motus, 
#  (iii) the corresponding values.

check2 <- reshape2::melt(NC_MiFish_clean$pcrs[,c("control_type", "nb_reads", "nb_motus")])

ggplot(data <- check2, aes(x=control_type, y=value, color=control_type)) + 
  geom_boxplot() + theme_bw() + 
  geom_jitter(alpha=0.2) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable, scales = "free_y") + 
  theme(axis.text.x = element_text(angle=45, h=1))

# Using the nb_reads and nb_motus defined previously in the NC_MiFish_clean$pcrs table

ggplot(NC_MiFish_clean$pcrs, aes(x=nb_reads, y=nb_motus, color = control_type)) + 
  geom_point() + theme_bw() + 
  scale_y_log10() + scale_x_log10() + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey")

ggpcrplate(NC_MiFish_clean, FUN = function(m){rowSums(m$reads)}, legend_title = "# of reads per PCR")

#check if this changes measures of beta diversity
# Get row data only for samples
tmp <- subset_metabarlist(NC_MiFish, table = "pcrs",
                          indices = NC_MiFish$pcrs$type == "sample" &
                            NC_MiFish$pcrs$nb_reads>0)

# Add sample biological information for checks
tmp$pcrs$material
#already exists so I don't need to grab a field from Samples for this 
tmp$pcrs$date_collected <- tmp$samples$`Date Collected`[match(tmp$pcrs$sample_id, rownames(tmp$samples))]
tmp$pcrs$date_collected
tmp$pcrs$species <- tmp$samples$`Species (common name)`[match(tmp$pcrs$sample_id, rownames(tmp$samples))]
tmp$pcrs$species

#NC_MiFish_clean$pcrs$material #already exists
NC_MiFish_clean$pcrs$date_collected <- NC_MiFish_clean$samples$`Date Collected`[match(NC_MiFish_clean$pcrs$sample_id,
                                                                                              rownames(NC_MiFish_clean$samples))]
NC_MiFish_clean$pcrs$species <-
  NC_MiFish_clean$samples$`Species (common name)`[match(NC_MiFish_clean$pcrs$sample_id,
                                                            rownames(NC_MiFish_clean$samples))]

#make tmp3 which is NC_MiFish_clean without any pcrs with 0 reads.
tmp3 <- subset_metabarlist(NC_MiFish_clean, table = "pcrs",
                           indices = #NC_MiFish_clean$pcrs$type == "sample" &
                             NC_MiFish_clean$pcrs$nb_reads>0)

# Add sample biological information for checks
tmp3$pcrs$material
#already exists so I don't need to grab a field from Samples for this 
tmp3$pcrs$date_collected <- tmp3$samples$`Date Collected`[match(tmp3$pcrs$sample_id, rownames(tmp3$samples))]
tmp3$pcrs$date_collected
tmp3$pcrs$species <- tmp3$samples$`Species (common name)`[match(tmp3$pcrs$sample_id, rownames(tmp3$samples))]
tmp3$pcrs$species


# Build PCoA ordinations 
mds1 <- check_pcr_repl(tmp,
                       groups = paste(
                         tmp$pcrs$species,
                         tmp$pcrs$material,
                         #tmp$pcrs$date_collected, 
                         sep = " | "))
mds2 <- check_pcr_repl(tmp3,
                       groups = paste(
                         tmp3$pcrs$species,
                         tmp3$pcrs$material,
                         #NC_MiFish_clean$pcrs$date_collected,
                         sep = " | "))


# Custom colors
a <- mds1 + labs(color = "species | material") + 
  scale_color_discrete()+#_manual(values = c("brown4", "brown1", "goldenrod4", "goldenrod1")) +
  theme(legend.position = "none") + 
  ggtitle("Raw data")
b <- mds2 + labs(color = "species | material") +
  scale_color_discrete()+#manual(values = c("brown4", "brown1", "goldenrod4", "goldenrod1")) + 
  ggtitle("Clean data")

# Assemble plots
leg <- get_legend(b + guides(shape=F) + 
                    theme(legend.position = "right", 
                          legend.direction = "vertical"))
ggdraw() +
  draw_plot(a, x=0, y=0, width = 0.4, height = 1) + 
  draw_plot(b + guides(color=F, shape=F), x=0.42, y=0, width = 0.4, height = 1) +
  draw_grob(leg, x=0.4, y=0)


#in the end I went with 0.04% cutoff for tagjumps and 98% similarity.

#write results to 4 csvs
write.csv(NC_MiFish_clean$reads, "data/06_MetabaR_results/NC_MiFish_obi3_lulu_metabar_reads.csv")
write.csv(NC_MiFish_clean$motus, "data/06_MetabaR_results/NC_MiFish_obi3_lulu_metabar_motus.csv")
write.csv(NC_MiFish_clean$pcrs, "data/06_MetabaR_results/NC_MiFish_obi3_lulu_metabar_pcrs.csv")
write.csv(NC_MiFish_clean$samples, "data/06_MetabaR_results/NC_MiFish_obi3_lulu_metabar_samples.csv")

#donor shark motus will make the different species fall out differently from each other on the PCA plot but removing them will help show the differences in their diets.



#remove the donor shark species from the results for each species
#donorspecies<-NC_MiFish_clean$samples$`Species (scientific name)`
#for (i in 1:length(NC_MiFish_clean$motus$species)) {
#  print(donorspecies)
#  NC_MiFish_clean$motus$donor_species[[i]] <- grepl(NC_MiFish_clean$motus$species[[i]], donorspecies)
#}

#NC_MiFish_clean$motus$donor_species <- 
#  ifelse(NC_MiFish_clean$motus$donor_species[match(NC_MiFish_clean$samples$`Species (scientific name)`, NC_MiFish_clean$motus$species)]!=TRUE,  F, T)

#NC_MiFish_clean$motus$donor_species <- sapply(NC_MiFish_clean$samples$`Species (scientific name)`,function(x) NC_MiFish_clean$motus$species[grep(x,NC_MiFish_clean$motus$species)[1]])

#df$V4 <- sapply(df$V3,function(x) df$V2[grep(x,df$V1)[1]])


# Proportion of each of these over total number of MOTUs
#table(NC_MiFish$motus$target_taxon) / nrow(NC_MiFish$motus)

# Intersection with extraction contaminant flags (not contaminant = T)
#table(NC_MiFish$motus$target_taxon, 
#      NC_MiFish$motus$not_a_pcr_conta)

NC_MiFish_clean$pcrs$species <- NC_MiFish_clean$samples$`Species (scientific name)`[match(NC_MiFish_clean$pcrs$sample_id, rownames(NC_MiFish_clean$samples))]

NC_MiFish_clean$motus$Sphyrna_tiburo_motus <- grepl("Sphyrna", NC_MiFish_clean$motus$SCIENTIFIC_NAME)
NC_MiFish_clean$motus$C_acronotus_motus <- grepl("Carcharhinus", NC_MiFish_clean$motus$SCIENTIFIC_NAME)
NC_MiFish_clean$motus$C_limbatus_motus <- grepl("Carcharhinus", NC_MiFish_clean$motus$SCIENTIFIC_NAME)
NC_MiFish_clean$motus$R_terranovae_motus <- grepl("Rhizoprionodon", NC_MiFish_clean$motus$SCIENTIFIC_NAME)

summary_metabarlist(NC_MiFish_clean)
#138 pcrs survived and 182 motus average 13 motus/sample


#A proxy for the individually separated ones, just to look at the PCA plots

NC_MiFish_no_donors <-subset_metabarlist(NC_MiFish_clean, table ="motus", indices =rowSums(NC_MiFish_clean$motus[,c("Sphyrna_tiburo_motus", "C_acronotus_motus","C_limbatus_motus", "R_terranovae_motus")]) == 0 )


# Compute the number of reads per pcr
NC_MiFish_no_donors$pcrs$nb_reads <- rowSums(NC_MiFish_no_donors$reads)

# Compute the number of motus per pcr
NC_MiFish_no_donors$pcrs$nb_motus <- rowSums(NC_MiFish_no_donors$reads>0)
summary_metabarlist(NC_MiFish_no_donors)
#136 motus 5.9 average motus/sample

NC_MiFish_no_donors1 <- subset_metabarlist(NC_MiFish_no_donors, table = "pcrs",
                           indices = NC_MiFish_no_donors$pcrs$nb_reads>0)



summary_metabarlist(NC_MiFish_no_donors1)

#when removing pcrs with 0 reads, we now have 81 samples surviving, and 136 motus average 9.9 motus/sample


# Build PCoA ordinations 
mds1 <- check_pcr_repl(tmp,
                       groups = paste(
                         tmp$pcrs$species,
                         tmp$pcrs$material,
                         sep = " | "))
mds4 <- check_pcr_repl(NC_MiFish_no_donors1,
                       groups = paste(
                         NC_MiFish_no_donors1$pcrs$species,
                         NC_MiFish_no_donors1$pcrs$material,
                         sep = " | "))



# Custom colors
a <- mds1 + labs(color = "species | material") + 
  scale_color_discrete()+#_manual(values = c("brown4", "brown1", "goldenrod4", "goldenrod1")) +
  theme(legend.position = "none") + 
  ggtitle("Raw data")
b <- mds4 + labs(color = "species | material") +
  scale_color_discrete()+#manual(values = c("brown4", "brown1", "goldenrod4", "goldenrod1")) + 
  ggtitle("Clean data")

# Assemble plots
leg <- get_legend(b + guides(shape=F) + 
                    theme(legend.position = "right", 
                          legend.direction = "vertical"))
ggdraw() +
  draw_plot(a, x=0, y=0, width = 0.4, height = 1) + 
  draw_plot(b + guides(color=F, shape=F), x=0.42, y=0, width = 0.4, height = 1) +
  draw_grob(leg, x=0.4, y=0)






#For subsetting the data
# Get the samples names from the stomachs plot
stomach_id <- grepl("St", rownames(NC_MiFish_clean$pcrs))

# Subset the data
NC_MiFish_St <- subset_metabarlist(NC_MiFish_clean, table = "pcrs", indices = stomach_id)

# Check results
summary_metabarlist(NC_MiFish_St)

#same for fecals
fecal_id<- grepl("NC[1-9]", rownames(NC_MiFish_clean$pcrs))
NC_MiFish_Fecals <- subset_metabarlist(NC_MiFish_clean, table = "pcrs", indices = fecal_id)  
summary_metabarlist(NC_MiFish_Fecals)

#now just subset the MiFish negative pcr controls
negative_control_id<- grepl("MF", rownames(NC_MiFish_clean$pcrs))
NC_MiFish_MF_Negatives <- subset_metabarlist(NC_MiFish_clean, table = "pcrs", indices = negative_control_id)  
summary_metabarlist(NC_MiFish_MF_Negatives)


#now just subset the negative pcr controls
negative_control_id<- grepl("neg", rownames(NC_MiFish_clean$pcrs))
NC_MiFish_Negatives <- subset_metabarlist(NC_MiFish_clean, table = "pcrs", indices = negative_control_id)  
summary_metabarlist(NC_MiFish_Negatives)

#View(NC_MiFish_Negatives[["motus"]])
#View(NC_MiFish_Negatives[["pcrs"]])
#View(NC_MiFish_Negatives[["reads"]])
#YAYAYAYAYAY zero motus or reads in the negatives.

#now make separate metabarlists for each species and remove just their own genus(only) or species DNA


#Blacknoses
C.acronotus_MiFish<- grepl("acronotus", NC_MiFish_clean$samples$`Species (scientific name)`)

NC_MiFish_C.acronotus <- subset_metabarlist(NC_MiFish_clean, table = "samples", indices = C.acronotus_MiFish)  
summary_metabarlist(NC_MiFish_C.acronotus)

NC_MiFish_C.acronotus_no_cannibals<- subset_metabarlist(NC_MiFish_C.acronotus, table ="motus", indices =rowSums(NC_MiFish_C.acronotus$motus["C_acronotus_motus"]) == 0 )
summary_metabarlist(NC_MiFish_C.acronotus_no_cannibals)

# Compute the number of reads per pcr
NC_MiFish_C.acronotus_no_cannibals$pcrs$nb_reads <- rowSums(NC_MiFish_C.acronotus_no_cannibals$reads)

# Compute the number of motus per pcr
NC_MiFish_C.acronotus_no_cannibals$pcrs$nb_motus <- rowSums(NC_MiFish_C.acronotus_no_cannibals$reads>0)
summary_metabarlist(NC_MiFish_C.acronotus_no_cannibals)
#110 motus 8 average motus/sample

blacknose <- subset_metabarlist(NC_MiFish_C.acronotus_no_cannibals, table = "pcrs",
                                           indices = NC_MiFish_C.acronotus_no_cannibals$pcrs$nb_reads>0)

summary_metabarlist(blacknose)


#110 motus, 21 sharks, 23 samples(pcrs in MetabaR terms), average 14 fish motus/sample

#Blacktips
C.limbatus_MiFish<- grepl("limbatus", NC_MiFish_clean$samples$`Species (scientific name)`)

NC_MiFish_C.limbatus <- subset_metabarlist(NC_MiFish_clean, table = "samples", indices = C.limbatus_MiFish)  
summary_metabarlist(NC_MiFish_C.limbatus)

NC_MiFish_C.limbatus_no_cannibals<- subset_metabarlist(NC_MiFish_C.limbatus, table ="motus", indices =rowSums(NC_MiFish_C.limbatus$motus["C_limbatus_motus"]) == 0 )
summary_metabarlist(NC_MiFish_C.limbatus_no_cannibals)

# Compute the number of reads per pcr
NC_MiFish_C.limbatus_no_cannibals$pcrs$nb_reads <- rowSums(NC_MiFish_C.limbatus_no_cannibals$reads)

# Compute the number of motus per pcr
NC_MiFish_C.limbatus_no_cannibals$pcrs$nb_motus <- rowSums(NC_MiFish_C.limbatus_no_cannibals$reads>0)
summary_metabarlist(NC_MiFish_C.limbatus_no_cannibals)
#80 motus 10 samples, 15 pcrs, 12.3 average motus/sample

blacktips <- subset_metabarlist(NC_MiFish_C.limbatus_no_cannibals, table = "pcrs",
                                indices = NC_MiFish_C.limbatus_no_cannibals$pcrs$nb_reads>0)

summary_metabarlist(blacktips)

#Atlantic sharpnose
R.terraenovae_MiFish<- grepl("Rhizoprionodon terraenovae", NC_MiFish_clean$samples$`Species (scientific name)`)

NC_MiFish_R.terraenovae <- subset_metabarlist(NC_MiFish_clean, table = "samples", indices = R.terraenovae_MiFish)  
summary_metabarlist(NC_MiFish_R.terraenovae)

NC_MiFish_R.terraenovae_no_cannibals<- subset_metabarlist(NC_MiFish_R.terraenovae, table ="motus", indices =rowSums(NC_MiFish_R.terraenovae$motus["R_terranovae_motus"]) == 0 )
summary_metabarlist(NC_MiFish_R.terraenovae_no_cannibals)

# Compute the number of reads per pcr
NC_MiFish_R.terraenovae_no_cannibals$pcrs$nb_reads <- rowSums(NC_MiFish_R.terraenovae_no_cannibals$reads)

# Compute the number of motus per pcr
NC_MiFish_R.terraenovae_no_cannibals$pcrs$nb_motus <- rowSums(NC_MiFish_R.terraenovae_no_cannibals$reads>0)
summary_metabarlist(NC_MiFish_R.terraenovae_no_cannibals)
#105 motus 39 samples, 45 pcrs, 5.5 average motus/sample

sharpnose <- subset_metabarlist(NC_MiFish_R.terraenovae_no_cannibals, table = "pcrs",
                                indices = NC_MiFish_R.terraenovae_no_cannibals$pcrs$nb_reads>0)

summary_metabarlist(sharpnose)


#105 motus 35 pcrs surviving (31 samples), 7.1 motus average/sample

#Bonnethead
S.tiburo_MiFish<- grepl("Sphyrna tiburo", NC_MiFish_clean$samples$`Species (scientific name)`)

NC_MiFish_S.tiburo <- subset_metabarlist(NC_MiFish_clean, table = "samples", indices = S.tiburo_MiFish)  
summary_metabarlist(NC_MiFish_S.tiburo)

NC_MiFish_S.tiburo_no_cannibals<- subset_metabarlist(NC_MiFish_S.tiburo, table ="motus", indices =rowSums(NC_MiFish_S.tiburo$motus["Sphyrna_tiburo_motus"]) == 0 )
summary_metabarlist(NC_MiFish_S.tiburo_no_cannibals)

# Compute the number of reads per pcr
NC_MiFish_S.tiburo_no_cannibals$pcrs$nb_reads <- rowSums(NC_MiFish_S.tiburo_no_cannibals$reads)

# Compute the number of motus per pcr
NC_MiFish_S.tiburo_no_cannibals$pcrs$nb_motus <- rowSums(NC_MiFish_S.tiburo_no_cannibals$reads>0)
summary_metabarlist(NC_MiFish_S.tiburo_no_cannibals)
#24 fish  motus 28 samples, 37 pcrs, 0.73 average fish motus/sample


bonnetheads <- subset_metabarlist(NC_MiFish_S.tiburo_no_cannibals, table = "pcrs",
                                indices = NC_MiFish_S.tiburo_no_cannibals$pcrs$nb_reads>0)

summary_metabarlist(bonnetheads)


#24 motus 11 pcrs surviving (10 samples), 2.45 fish motus average/sample

#Write 4 separate files for each donor species

#blacknose
write.csv(NC_MiFish_C.acronotus_no_cannibals$reads, "data/06_MetabaR_results/NC_MiFish_C.acronotus_no_cannibals_metabar_reads.csv")
write.csv(NC_MiFish_C.acronotus_no_cannibals$motus, "data/06_MetabaR_results/NC_MiFish_C.acronotus_no_cannibals_metabar_motus.csv")
write.csv(NC_MiFish_C.acronotus_no_cannibals$pcrs, "data/06_MetabaR_results/NC_MiFish_C.acronotus_no_cannibals_metabar_pcrs.csv")
write.csv(NC_MiFish_C.acronotus_no_cannibals$samples, "data/06_MetabaR_results/NC_MiFish_C.acronotus_no_cannibals_metabar_samples.csv")




#blacktips
write.csv(NC_MiFish_C.limbatus_no_cannibals$reads, "data/06_MetabaR_results/NC_MiFish_C.limbatus_no_cannibals_metabar_reads.csv")
write.csv(NC_MiFish_C.limbatus_no_cannibals$motus, "data/06_MetabaR_results/NC_MiFish_C.limbatus_no_cannibals_metabar_motus.csv")
write.csv(NC_MiFish_C.limbatus_no_cannibals$pcrs, "data/06_MetabaR_results/NC_MiFish_C.limbatus_no_cannibals_metabar_pcrs.csv")
write.csv(NC_MiFish_C.limbatus_no_cannibals$samples, "data/06_MetabaR_results/NC_MiFish_C.limbatus_no_cannibals_metabar_samples.csv")

#sharpnose
write.csv(NC_MiFish_R.terraenovae_no_cannibals$reads, "data/06_MetabaR_results/NC_MiFish_R.terraenovae_no_cannibals_metabar_reads.csv")
write.csv(NC_MiFish_R.terraenovae_no_cannibals$motus, "data/06_MetabaR_results/NC_MiFish_R.terraenovae_no_cannibals_metabar_motus.csv")
write.csv(NC_MiFish_R.terraenovae_no_cannibals$pcrs, "data/06_MetabaR_results/NC_MiFish_R.terraenovae_no_cannibals_metabar_pcrs.csv")
write.csv(NC_MiFish_R.terraenovae_no_cannibals$samples, "data/06_MetabaR_results/NC_MiFish_R.terraenovae_no_cannibals_metabar_samples.csv")


#bonnetheads
write.csv(NC_MiFish_S.tiburo_no_cannibals$reads, "data/06_MetabaR_results/NC_MiFish_S.tiburo_no_cannibals_metabar_reads.csv")
write.csv(NC_MiFish_S.tiburo_no_cannibals$motus, "data/06_MetabaR_results/NC_MiFish_S.tiburo_no_cannibals_metabar_motus.csv")
write.csv(NC_MiFish_S.tiburo_no_cannibals$pcrs, "data/06_MetabaR_results/NC_MiFish_S.tiburo_no_cannibals_metabar_pcrs.csv")
write.csv(NC_MiFish_S.tiburo_no_cannibals$samples, "data/06_MetabaR_results/NC_MiFish_S.tiburo_no_cannibals_metabar_samples.csv")

#combine all MiFish "no cannibals" metabar reads files with the final BerryCrust reads file

df.S.tiburo_no_cannibals_reads<- rownames_to_column(as.data.frame(NC_MiFish_S.tiburo_no_cannibals$reads), var="sample")
df.R.terraenovae_no_cannibals_reads<- rownames_to_column(as.data.frame(NC_MiFish_R.terraenovae_no_cannibals$reads), var="sample")
df.C.limbatus_no_cannibals_reads<- rownames_to_column(as.data.frame(NC_MiFish_C.limbatus_no_cannibals$reads), var="sample")
df.C.acronotus_no_cannibals_reads<- rownames_to_column(as.data.frame(NC_MiFish_C.acronotus_no_cannibals$reads), var="sample")

#combine all no cannibals reads

df.MiFish.no_cannibals.reads<-full_join(full_join(full_join(df.S.tiburo_no_cannibals_reads,df.R.terraenovae_no_cannibals_reads), df.C.limbatus_no_cannibals_reads),df.C.acronotus_no_cannibals_reads)

#get the results from the BerryCrust NC_Sharks R script with Metabar
df.BerryCrust_results.reads<-read.csv("data/06_MetabaR_results/NC_BerryCrust_obi3_lulu_metabar_reads.csv")
df.BerryCrust_results.reads<-rename(df.BerryCrust_results.reads, "X"="sample")

df.NC_Sharks_results_merged_reads<-full_join(df.MiFish.no_cannibals.reads,df.BerryCrust_results.reads)

df.NC_Sharks_results_merged_reads1<-column_to_rownames(df.NC_Sharks_results_merged_reads, var="sample")
write.csv(df.NC_Sharks_results_merged_reads1, "NC_Sharks_merged_final_reads.csv")

#motus
df.S.tiburo_no_cannibals_motus<- rownames_to_column(as.data.frame(NC_MiFish_S.tiburo_no_cannibals$motus), var="id")
df.R.terraenovae_no_cannibals_motus<- rownames_to_column(as.data.frame(NC_MiFish_R.terraenovae_no_cannibals$motus), var="id")
df.C.limbatus_no_cannibals_motus<- rownames_to_column(as.data.frame(NC_MiFish_C.limbatus_no_cannibals$motus), var="id")
df.C.acronotus_no_cannibals_motus<- rownames_to_column(as.data.frame(NC_MiFish_C.acronotus_no_cannibals$motus), var="id")

#combine all no cannibals motus

df.MiFish.no_cannibals.motus<-full_join(full_join(full_join(df.S.tiburo_no_cannibals_motus,df.R.terraenovae_no_cannibals_motus), df.C.limbatus_no_cannibals_motus),df.C.acronotus_no_cannibals_motus)

#get the results from the BerryCrust NC_Sharks R script with Metabar
df.BerryCrust_results.motus<-read.csv("data/06_MetabaR_results/NC_BerryCrust_obi3_lulu_metabar_motus.csv")
df.BerryCrust_results.motus<-rename(df.BerryCrust_results.motus, "X"="id")

df.NC_Sharks_results_merged_motus<-full_join(df.MiFish.no_cannibals.motus,df.BerryCrust_results.motus)

df.NC_Sharks_results_merged_motus1<-column_to_rownames(df.NC_Sharks_results_merged_motus, var="id")
write.csv(df.NC_Sharks_results_merged_motus1, "NC_Sharks_merged_final_motus.csv")


#pcrs
df.S.tiburo_no_cannibals_pcrs<- rownames_to_column(as.data.frame(NC_MiFish_S.tiburo_no_cannibals$pcrs), var="sample")
df.R.terraenovae_no_cannibals_pcrs<- rownames_to_column(as.data.frame(NC_MiFish_R.terraenovae_no_cannibals$pcrs), var="sample")
df.C.limbatus_no_cannibals_pcrs<- rownames_to_column(as.data.frame(NC_MiFish_C.limbatus_no_cannibals$pcrs), var="sample")
df.C.acronotus_no_cannibals_pcrs<- rownames_to_column(as.data.frame(NC_MiFish_C.acronotus_no_cannibals$pcrs), var="sample")

#combine all no cannibals pcrs

df.MiFish.no_cannibals.pcrs<-full_join(full_join(full_join(df.S.tiburo_no_cannibals_pcrs,df.R.terraenovae_no_cannibals_pcrs), df.C.limbatus_no_cannibals_pcrs),df.C.acronotus_no_cannibals_pcrs)

#get the results from the BerryCrust NC_Sharks R script with Metabar
df.BerryCrust_results.pcrs<-read.csv("data/06_MetabaR_results/NC_BerryCrust_obi3_lulu_metabar_pcrs.csv")
df.BerryCrust_results.pcrs<-rename(df.BerryCrust_results.pcrs, "X"="sample")

df.NC_Sharks_results_merged_pcrs<-full_join(df.MiFish.no_cannibals.pcrs,df.BerryCrust_results, by=c("sample","plate_no","plate_row","tag_fwd","tag_rev","project","sample_id","material","type","control_type"))

df.NC_Sharks_results_merged_pcrs1<-column_to_rownames(df.NC_Sharks_results_merged_pcrs, var="sample")
write.csv(df.NC_Sharks_results_merged_pcrs1, "NC_Sharks_merged_final_pcrs.csv")



#samples
df.S.tiburo_no_cannibals_samples<- rownames_to_column(as.data.frame(NC_MiFish_S.tiburo_no_cannibals$samples), var="sample_id")
df.R.terraenovae_no_cannibals_samples<- rownames_to_column(as.data.frame(NC_MiFish_R.terraenovae_no_cannibals$samples), var="sample_id")
df.C.limbatus_no_cannibals_samples<- rownames_to_column(as.data.frame(NC_MiFish_C.limbatus_no_cannibals$samples), var="sample_id")
df.C.acronotus_no_cannibals_samples<- rownames_to_column(as.data.frame(NC_MiFish_C.acronotus_no_cannibals$samples), var="sample_id")

#combine all no cannibals reads

df.MiFish.no_cannibals.samples<-full_join(full_join(full_join(df.S.tiburo_no_cannibals_samples,df.R.terraenovae_no_cannibals_samples), df.C.limbatus_no_cannibals_samples),df.C.acronotus_no_cannibals_samples)

#get the results from the BerryCrust NC_Sharks R script with Metabar
df.BerryCrust_results.samples<-read_csv("data/06_MetabaR_results/NC_BerryCrust_obi3_lulu_metabar_samples.csv")
df.BerryCrust_results.samples<-rename(df.BerryCrust_results.samples, "...1"="sample_id")

df.BerryCrust_results.samples$`Time of Capture`<- as.character(df.BerryCrust_results.samples$`Time of Capture`)

df.BerryCrust_results.samples$`Time placed in incubator for Step 2`<- as.character(df.BerryCrust_results.samples$`Time placed in incubator for Step 2`)


df.NC_Sharks_results_merged_samples<-union_all(df.MiFish.no_cannibals.samples,df.BerryCrust_results.samples)
df.NC_Sharks_results_merged_samples1<-df.NC_Sharks_results_merged_samples %>% distinct(sample_id, .keep_all = TRUE)


df.NC_Sharks_results_merged_samples1<-column_to_rownames(df.NC_Sharks_results_merged_samples1, var="sample_id")
write.csv(df.NC_Sharks_results_merged_samples1, "data/06_MetabaR_results/NC_Sharks_merged_final_samples.csv")



#bring back into metabar for the full dataset PCA plot.
NC_Combined<-tabfiles_to_metabarlist(
  file_reads= "data/06_MetabaR_results/NC_Sharks_merged_final_reads_NAto0.csv",
  file_motus= "data/06_MetabaR_results/NC_Sharks_merged_final_motus_cleared_tags.csv",
  file_pcrs ="data/06_MetabaR_results/NC_Sharks_merged_final_pcrs.csv",
  file_samples ="data/06_MetabaR_results/NC_Sharks_merged_final_samples.csv",
  files_sep = ","
)


#Error in check_metabarlist(out) : 
#columns of table `reads` in out must correspond exactly to rows of `motus`

#it's the same number of each so it must be the sorting

#opened these files in excel and looked at the order, they look like they're ordered the same but I see NAs instead of 0s. replaced those. saved as ...reads _NAto0.csv and motus_cleared_tags to put the sequence last.
#!Worked! wth warnings about PCR tags and #Some PCRs in out have a number of reads of zero in table `reads`!

#PCA PLOT OF COMPLETE COMBINED DATASET!


# Compute the number of reads per pcr
NC_Combined$pcrs$nb_reads <- rowSums(NC_Combined$reads)

# Compute the number of motus per pcr
NC_Combined$pcrs$nb_motus <- rowSums(NC_Combined$reads>0)

df.NC_Combined_reads<-as.data.frame(NC_Combined$reads)


summary_metabarlist(NC_Combined)
#283 motus  121 sharks 146 pcr samples, 12.9 average motus/sample

NC_Combined1 <- subset_metabarlist(NC_Combined, table = "pcrs",
                                           indices = NC_Combined$pcrs$nb_reads>0)

summary_metabarlist(NC_Combined1)

#when removing pcrs with 0 reads, we now have 131 pcr samples surviving from 108 sharks,and 283 motus. average 14.4motus/sample


# Build PCoA ordination

NC_Combined1$pcrs$species<-NC_Combined1$samples$`Species (common name)`[match(NC_Combined1$pcrs$sample_id,rownames(NC_Combined1$samples))]

mds5 <- check_pcr_repl(NC_Combined1,
                       groups = paste(
                         NC_Combined1$pcrs$species,
                         NC_Combined1$pcrs$material,
                         sep = " | "))



# Custom colors
c <- mds5 + labs(color = "species | material") + 
  scale_color_discrete()+#_manual(values = c("brown4", "brown1", "goldenrod4", "goldenrod1")) +
  theme(legend.position = "none") + 
  ggtitle("Differences in Shark Diets by Species and Material")

# Assemble plots
leg <- get_legend(c + guides(shape=F) + 
                    theme(legend.position = "right", 
                          legend.direction = "vertical"))
ggdraw() +
  draw_plot(c + guides(color=F, shape=F), x=0, y=0, width = .75, height = 1) +
  draw_grob(leg, x=.4, y=0)

ggsave("data/07_Phyloseq_results_and_Visualizations/PCA_of_NC_Sharks_Diets.jpg")



#Save R Global Environment
#save.image(file = "data/06_MetabaR_results/Combined_MiFish_and_NC_Sharks_Rglobalenvironment.RData")





#Full Dataset Phyloseq starting from existing R objects made by the above code:


#Phyloseq!
#otu_table - Works on any numeric matrix. You must also specify if the species are rows or columns
#sample_data - Works on any data.frame. The rownames must match the sample names in the otu_table if you plan to combine them as a phyloseq-object
#tax_table - Works on any character matrix. The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
#phyloseq - Takes as argument an otu_table and any unordered list of valid phyloseq components: sample_data, tax_table, phylo, or XStringSet. The tip labels of a phylo-object (tree) must match the OTU names of the otu_table, and similarly, the sequence names of an XStringSet object must match the OTU names of the otu_table.
#merge_phyloseq - Can take any number of phyloseq objects and/or phyloseq components, and attempts to combine them into one larger phyloseq object. This is most-useful for adding separately-imported components to an already-created phyloseq object.


load("data/06_MetabaR_results/Combined_MiFish_and_NC_Sharks_Rglobalenvironment.RData")
#otu table

clean.otu.df<-as.data.frame(t(NC_Combined1$reads))
clean.otu.df<-rownames_to_column(clean.otu.df, var="id")

#taxa table

library(taxonomizr)
prepareDatabase(getAccessions=FALSE)
      #make a taxa table for all obi3 result otus
taxaId<-df.NC_Sharks_results_merged_motus$TAXID
taxa<-getTaxonomy(taxaId,'nameNode.sqlite')
print(taxa)
class(taxa)
taxa.df<-as.data.frame.array(taxa)
rownames(taxa.df)<-df.NC_Sharks_results_merged_motus$id

#samples table
df.NC_Sharks_results_merged_samples<-rownames_to_column(df.NC_Sharks_results_merged_samples1, var="sample_id")

samples.df<-full_join(df.NC_Sharks_results_merged_pcrs,df.NC_Sharks_results_merged_samples, by="sample_id")


samples.df <- samples.df %>% 
  tibble::column_to_rownames("sample") 

#samples.df<-samples.df[c(1,10,32:61)]

#make otu matrix
clean.otu.df<-column_to_rownames(clean.otu.df, var="id")
clean.otu.mat<-as.matrix(clean.otu.df)

#make taxa matrix
clean.taxa.mat<-as.matrix(taxa.df)

#Save and Load R Global Environment
#save.image(file = "data/06_MetabaR_results/NC_Sharks_ready_for_phyloseq.RData")
load("data/06_MetabaR_results/NC_Sharks_ready_for_phyloseq.RData")

library(phyloseq)
#Starting Phlyoseq!
NC_Sharks.ps<- phyloseq(otu_table(clean.otu.mat, taxa_are_rows = TRUE), 
                         sample_data(samples.df), 
                         tax_table(clean.taxa.mat))

NC_Sharks.ps
sample_names(NC_Sharks.ps)
rank_names(NC_Sharks.ps)
sample_variables(NC_Sharks.ps)

write.csv(clean.otu.df, "data/06_MetabaR_results/combined_NC_Sharks_ready_for_phyloseq_otus_no_empty_pcrs.csv")
write.csv(samples.df, "data/06_MetabaR_results/combined_NC_Sharks_ready_for_phyloseq_samples.csv")
write.csv(taxa.df, "data/06_MetabaR_results/combined_NC_Sharks_ready_for_phyloseq_taxa.csv")


NC_Sharks.ps
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 287 taxa and 131 samples ]
#sample_data() Sample Data:       [ 131 samples by 61 sample variables ]
#tax_table()   Taxonomy Table:    [ 287 taxa by 7 taxonomic ranks ]
#when it says taxa it means motus

plot_richness(NC_Sharks.ps, x="material", measures=c("Shannon", "Simpson"), color="material")

# Transform data to proportions as appropriate for Bray-Curtis distances
NC_Sharks.ps.prop <- transform_sample_counts(NC_Sharks.ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(NC_Sharks.ps.prop, method="NMDS", distance="bray")

#Warning message:
#In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data

plot_ordination(NC_Sharks.ps.prop, ord.nmds.bray, color="Species..scientific.name."
                , title="Bray NMDS")
#what does this mean?! changed color to species.

top20 <- names(sort(taxa_sums(NC_Sharks.ps), decreasing=TRUE))[0:20]
ps.top20 <- transform_sample_counts(NC_Sharks.ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Species..scientific.name.",fill="genus")
ggsave("Bar_Chart_of_Species_by_top20preyGenus.jpg")
###############

#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(NC_Sharks.ps))
standf = function(x, t=total) round(t * (x / sum(x)))
NC_Sharks.ps = transform_sample_counts(NC_Sharks.ps, standf)

plot_bar(NC_Sharks.ps, x="Species..scientific.name.",fill="genus")
ggsave("Bar_Chart_of_Species_by_PreyGenus.jpg")

plot_bar(NC_Sharks.ps, fill = "genus", facet_grid = "Species..scientific.name.")
ggsave("Faceted_sample_species_prey-genus.jpg")

plot_bar(NC_Sharks.ps, fill = "family", facet_grid = "Species..scientific.name.")
ggsave("Faceted_sample_species_prey-family.jpg")
