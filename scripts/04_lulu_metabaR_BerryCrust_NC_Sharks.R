#phyloseq of lulu curated obitools3 results 
#12-21-22

#install.packages("devtools")
library(devtools)
#install_github("tobiasgf/lulu")
require(lulu)
library(phyloseq)

setwd("~/Desktop/Savannah/following_obitools3_method/results_obi3-12-21-22")


#First, lulu curation of obi3 results file.

#(not matching merged obi3 sequences back to the results for the lulu step, just results-results for matchlist and results otutable)
otutab <- read.csv("BerryCrust_named_tab_LULU.txt", sep='\t', header=TRUE, as.is=TRUE, row.names = 1) #was made using the obi3 merged reads matched against the obi3 results 
matchlist <- read.table("BerryCrust_named_matchlist.txt", header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)
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
write.csv(curated_result$curated_table,"lulu_curated_BerryCrust_named.tab")


#Phyloseq section
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
#import data:

#otus
BerryCrust_lulu_table<-otu_table(curated_result$curated_table, taxa_are_rows = TRUE)
Metabar_formatted_samples<-read.csv("../MetabaR/MetabaR_formatted_samples_NCSharks.txt", sep="\t")

#samples
Metabar_formatted_pcrs<-read.csv("../MetabaR/MetabaR_formatted_pcrs_NCSharks_BerryCrust.txt", sep="\t")
samples.df<-full_join(Metabar_formatted_pcrs,Metabar_formatted_samples,by="sample_id")

BerryCrust_sample_data<-sample_data(samples.df)

#taxonomy
BerryCrust_named_obi3_results<-read.csv("BerryCrust_named.tab", sep="\t")
BerryCrust_named_obi3_results<-BerryCrust_named_obi3_results%>% select(id, BEST_IDENTITY, TAXID, SCIENTIFIC_NAME, ID_STATUS)

#install.packages("taxonomizr")
library(taxonomizr)
prepareDatabase(getAccessions=FALSE)
#make a taxa table for all obi3 result otus
taxaId<-BerryCrust_named_obi3_results$TAXID
taxa<-getTaxonomy(taxaId,'nameNode.sqlite')
print(taxa)
class(taxa)
taxa.df<-as.data.frame.array(taxa)
rownames(taxa.df)<-BerryCrust_named_obi3_results$id

#make a taxa table for only the lulu-curated otus
BerryCrust_lulu_table.df<-as.data.frame(BerryCrust_lulu_table)
BerryCrust_lulu_table.df<-rownames_to_column(BerryCrust_lulu_table.df, var ="id")
taxa.df<-rownames_to_column(taxa.df, var="id")
taxa.df.lulu<- semi_join(taxa.df,BerryCrust_lulu_table.df,by="id")

#see which taxa were discarded by lulu
#taxa.df.discards<-anti_join(taxa.df,BerryCrust_lulu_table.df,by="id")

#do something about the negatives

#Metabar!

#format the motus
#taxa.df.lulu has the right motus and taxonomic levels, BerryCrust_named_obi3_motus has the rest of the info I want but too many motus join by ="id"

Metabar_formatted_motus<- inner_join(BerryCrust_named_obi3_results,taxa.df.lulu,by="id")
Metabar_formatted_motus<-column_to_rownames(Metabar_formatted_motus, var="id")


#flip rows and columns (transpose) of BerryCrust_lulu_table.df to make the reads table
BerryCrust_lulu_table.df<-as.data.frame(BerryCrust_lulu_table)

Metabar_formatted_reads<-t(BerryCrust_lulu_table.df)

Metabar_formatted_reads<-as.data.frame(Metabar_formatted_reads)

Metabar_formatted_reads<-as.matrix(Metabar_formatted_reads)

Metabar_formatted_samples<-column_to_rownames(Metabar_formatted_samples, var="sample_id")

Metabar_formatted_pcrs<-column_to_rownames(Metabar_formatted_pcrs, var="sample")



#NC_BerryCrust<-metabarlist_generator(reads=Metabar_formatted_reads, Metabar_formatted_motus, pcrs=Metabar_formatted_pcrs, samples=Metabar_formatted_samples)
#this seems to work, but doesn't identify any pcr contaminants which makes me distrust it, so, I'm going to export to csv and look at them in excel and reimport with tabfilestometabarlist

#write.csv(Metabar_formatted_reads, "Metabar_formatted_reads_pureobi3_lulu.csv")
#write.csv(Metabar_formatted_motus, "Metabar_formatted_motus_pureobi3_lulu.csv")
#write.csv(Metabar_formatted_pcrs, "Metabar_formatted_pcrs_pureobi3_lulu.csv")
#write.csv(Metabar_formatted_samples, "Metabar_formatted_samples_pureobi3_lulu.csv")


#samples file has lots of extra periods in it and the first line doesn't say sample_id, so I saved MetabaR_formatted_samples_NCSharks.txt to this folder to use instead.
#pcrs file looked exactly the same as a previous one I used for MiFish analysis except the first column label was gone so I put it back and saved it as a tab delimited text. Metabar_formatted_pcrs_pureobi3_lulu.txt
#motus needed some rearranging to make it like one of the MiFish ones that worked in the past.  I put the count field second and made it "count" instead of "COUNT", then the second field is seq_length and I calculated that in excel from the sequence field, which I put at the end of the taxonomy fields. Saved as Metabar_formatted_motus_pureobi3_lulu.txt
#reads all looked good so I saved as Metabar_formatted_reads_pureobi3_lulu.txt

#tabfiles_to_metabarlist(
#  file_reads= "Metabar_formatted_reads_pureobi3_lulu.txt",
#  file_motus= "Metabar_formatted_motus_pureobi3_lulu.txt",
#  file_pcrs ="Metabar_formatted_pcrs_pureobi3_lulu.txt",
#  file_samples ="MetabaR_formatted_samples_NCSharks.txt",
#  files_sep = "\t"
#)
``#Error in tabfiles_to_metabarlist(file_reads = "Metabar_formatted_reads_pureobi3_lulu.txt",  : 
#cannot continue, rownames in reads are not part of rownames of pcrs

#looked at these two and added sample to the beginning of the reads file to match PCRs, and sorted them both alphabetically so they match.  Motus wasn't sorted alphabetically to begin with.

#same error but pcrs had St9 and reads did not.  So I removed the St9 row from pcrs and saved as Metabar_formatted_pcrs_pureobi3_lulu_noSt9.txt

NC_BerryCrust<-tabfiles_to_metabarlist(
  file_reads= "Metabar_formatted_reads_pureobi3_lulu.txt",
  file_motus= "Metabar_formatted_motus_pureobi3_lulu.txt",
  file_pcrs ="Metabar_formatted_pcrs_pureobi3_lulu_noSt9.txt",
  file_samples ="MetabaR_formatted_samples_NCSharks.txt",
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
NC_BerryCrust$pcrs$nb_reads <- rowSums(NC_BerryCrust$reads)

# Compute the number of motus per pcr
NC_BerryCrust$pcrs$nb_motus <- rowSums(NC_BerryCrust$reads>0)

# Load requested package for plotting
library(ggplot2)
library(reshape2)
library(tidyverse)

# Create an input table (named check1) for ggplot of 3 columns: 
#  (i) control type 
#  (ii) a vector indicated whether it corresponds to nb_reads or nb_motus, 
#  (iii) the corresponding values.

check1 <- reshape2::melt(NC_BerryCrust$pcrs[,c("control_type", "nb_reads", "nb_motus")])

ggplot(data <- check1, aes(x=control_type, y=value, color=control_type)) + 
  geom_boxplot() + theme_bw() + 
  geom_jitter(alpha=0.2) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable, scales = "free_y") + 
  theme(axis.text.x = element_text(angle=45, h=1))

# Using the nb_reads and nb_motus defined previously in the NC_BerryCrust$pcrs table

ggplot(NC_BerryCrust$pcrs, aes(x=nb_reads, y=nb_motus, color = control_type)) + 
  geom_point() + theme_bw() + 
  scale_y_log10() + scale_x_log10() + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey")

ggpcrplate(NC_BerryCrust, FUN = function(m){rowSums(m$reads)}, legend_title = "# of reads per PCR")


# Here the list of all tag/indices used in the experiment 
# is available in the column "tag_rev" of the NC_MiFish$pcrs table
tag.list <- as.character(unique(NC_BerryCrust$pcrs$tag_fwd))
ggpcrtag(NC_BerryCrust, 
         legend_title = "# of reads per PCR", 
         FUN = function(m) {rowSums(m$reads)},
         taglist = tag.list) 


NC_BerryCrust.raref = hill_rarefaction(NC_BerryCrust, nboot = 20, nsteps = 10)
head(NC_BerryCrust.raref$hill_table)
gghill_rarefaction(NC_BerryCrust.raref) 


# Define a vector containing the Material info for each pcrs 
material <- NC_BerryCrust$pcrs$material

# Use of gghill_rarefaction requires a vector with named pcrs
material <- setNames(material,rownames(NC_BerryCrust$pcrs))

# Plot
p <- gghill_rarefaction(NC_BerryCrust.raref, group=material)
p + scale_fill_manual(values = c("goldenrod4", "brown4", "grey")) +
  scale_color_manual(values = c("goldenrod4", "brown4", "grey")) +
  labs(color="Material type")


# Identifying extraction contaminants(changed to pcr since that's the only control type we have)
NC_BerryCrust <- contaslayer(NC_BerryCrust, 
                         control_types = "pcr",
                         output_col = "not_a_pcr_conta")

table(NC_BerryCrust$motus$not_a_pcr_conta)
#I blasted the top 3 contaminants and they came back with very low confidence in catsharks and zebrasharks.

# Identify the most common contaminant
# get contaminant ids
id <- !NC_BerryCrust$motus$not_a_pcr_conta
max.conta <- rownames(NC_BerryCrust$motus[id,])[which.max(NC_BerryCrust$motus[id, "count"])]

#... and its distribution and relative abundance in each pcr
ggpcrplate(NC_BerryCrust, legend_title = "#reads of most \nabundant contaminant",
           FUN = function(m) {m$reads[, max.conta]/rowSums(m$reads)})

# Compute relative abundance of all pcr contaminants together 
a <- data.frame(conta.relab = rowSums(NC_BerryCrust$reads[,!NC_BerryCrust$motus$not_a_pcr_conta]) / 
                  rowSums(NC_BerryCrust$reads))
# Add information on control types
a$control_type <- NC_BerryCrust$pcrs$control_type[match(rownames(a), rownames(NC_BerryCrust$pcrs))]

ggplot(a, aes(x=control_type, y=conta.relab, color=control_type)) + 
  geom_boxplot() + geom_jitter(alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x=NULL, y="Prop. Reads (log10)") + 
  theme_bw() + 
  scale_y_log10()


# flag pcrs with total contaminant relative abundance > 10% of reads)
NC_BerryCrust$pcrs$low_contamination_level <- 
  ifelse(a$conta.relab[match(rownames(NC_BerryCrust$pcrs), rownames(a))]>1e-1,  F, T)

# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(NC_BerryCrust$pcrs$low_contamination_level) / nrow(NC_BerryCrust$pcrs)

#next the tutorial excludes non-eukaryotes and plots the similarity scores and removes the ones with low scores.  
#Whatever metabarcoding pipeline hasn't done that yet isn't worth it's salt IMO

#Flag MOTUs corresponding to target (TRUE) vs. non-target (FALSE) taxa 
NC_BerryCrust$motus$target_taxon <- grepl("TRUE", NC_BerryCrust$motus$ID_STATUS)#I changed this from Eukaryota and path because I am using the obitools3 results instead of the sintax vsearch results

# Proportion of each of these over total number of MOTUs
table(NC_BerryCrust$motus$target_taxon) / nrow(NC_BerryCrust$motus)

# Intersection with extraction contaminant flags (not contaminant = T)
table(NC_BerryCrust$motus$target_taxon, 
      NC_BerryCrust$motus$not_a_pcr_conta)


#remove non-target taxa and BC control contaminants
NC_BerryCrust <- contaslayer(NC_BerryCrust, 
                             method= "all",
                             control_types = "pcr",
                             controls = c("NCFSBCneg","NCSTBCneg"),
                             output_col = "not_a_BC_pcr_conta")

#table(NC_BerryCrust$motus$not_a_pcr_conta)
#I'd like to be able to look at which 9 were FALSE
table(NC_BerryCrust$motus$not_a_BC_pcr_conta)
table(NC_BerryCrust$motus$target_taxon, 
      NC_BerryCrust$motus$not_a_BC_pcr_conta)
table(NC_BerryCrust$motus$not_a_BC_pcr_conta, 
      NC_BerryCrust$motus$not_a_pcr_conta)

#all of the pcr contas are BC_pcr_contas, and none are target taxa! yay!

NC_BerryCrust1 <- subset_metabarlist(NC_BerryCrust, "motus", 
                          indices = rowSums(NC_BerryCrust$motus[,c("not_a_BC_pcr_conta", "target_taxon")]) == 2)


summary_metabarlist(NC_BerryCrust)

summary_metabarlist(NC_BerryCrust1)


# Plot the unweighted distribution of MOTUs similarity scores 
a <- 
  ggplot(NC_BerryCrust1$motus, aes(x=BEST_IDENTITY)) + 
  geom_histogram(color="grey", fill="white", bins=20) + 
  geom_vline(xintercept = 0.97, col="orange", lty=2) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(x="% similarity against best match", y="# MOTUs")

# Same for the weighted distribution
b <- 
  ggplot(NC_BerryCrust1$motus, 
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
NC_BerryCrust1$motus$not_degraded <-
  ifelse(NC_BerryCrust1$motus$BEST_IDENTITY < 0.98, F, T)

# Proportion of each of these over total number of MOTUs
table(NC_BerryCrust1$motus$not_degraded) / nrow(NC_BerryCrust1$motus)

# Flag not degraded (TRUE) vs. potentially degraded sequences (FALSE)
NC_BerryCrust$motus$not_degraded <-
  ifelse(NC_BerryCrust$motus$BEST_IDENTITY < 0.98, F, T)

# Proportion of each of these over total number of MOTUs
table(NC_BerryCrust$motus$not_degraded) / nrow(NC_BerryCrust$motus)
# for NC_BerryCrust
#> 
#    FALSE      TRUE 
#0.8042813 0.1957187 
#> 
#> For NC_BerryCrust1
#
#FALSE      TRUE 
#0.1351351 0.8648649 

# Intersection with other flags
table(NC_BerryCrust$motus$target_taxon, 
      NC_BerryCrust$motus$not_a_pcr_conta, 
      NC_BerryCrust$motus$not_degraded)


#THIS LEAVES 20 MOTUS! CHECK THE TAXA AND SEE IF A 97% SIMILARITY VALUE IS BETTER.

#Moving on to detecting PCR outliers by sequencing depth
ggplot(NC_BerryCrust$pcrs, aes(nb_reads)) +
  geom_histogram(bins=40, color="grey", fill="white") + 
  geom_vline(xintercept = 2e3, lty=2, color="orange") + # threshold
  scale_x_log10() + 
  labs(x="# Reads (with all MOTUs and PCRs)", 
       y="# PCRs") +
  theme_bw() + 
  theme(panel.grid = element_blank())


#play with the xintercept above to visualize your cutoff level

# Flag pcrs with an acceptable sequencing depth (TRUE) or inacceptable one (FALSE)
NC_BerryCrust$pcrs$seqdepth_ok <- ifelse(NC_BerryCrust$pcrs$nb_reads < 2e3, F, T)

# Proportion of each of these over total number of pcrs, control excluded
table(NC_BerryCrust$pcrs$seqdepth_ok[NC_BerryCrust$pcrs$type=="sample"]) /
  nrow(NC_BerryCrust$pcrs[NC_BerryCrust$pcrs$type=="sample",])
#and for the semi-cleaned dataset
NC_BerryCrust1$pcrs$seqdepth_ok <- ifelse(NC_BerryCrust1$pcrs$nb_reads < 2e3, F, T)

# Proportion of each of these over total number of pcrs, control excluded
table(NC_BerryCrust1$pcrs$seqdepth_ok[NC_BerryCrust1$pcrs$type=="sample"]) /
  nrow(NC_BerryCrust1$pcrs[NC_BerryCrust1$pcrs$type=="sample",])


#They use this next function to compare the similarity of biological controls with the expectation they will be more similar to each other than to other samples.  
#We can use this to check the similarity of fecals to their associated stomachs

# Subsetting the metabarlist
NC_BerryCrust_sub <- subset_metabarlist(NC_BerryCrust, 
                                    table = "pcrs", 
                                    indices = NC_BerryCrust$pcrs$nb_reads>0 & (
                                      is.na(NC_BerryCrust$pcrs$control_type) |
                                        NC_BerryCrust$pcrs$control_type=="positive"))

# First visualization
comp1 = pcr_within_between(NC_BerryCrust_sub)
check_pcr_thresh(comp1)

#cool, this looks good.  Now... for flagging not-well-replicated pcrs.  Probably shouldn't flag these for removal since that's not appropriate for this experimental design

# Subsetting the metabarlist
NC_BerryCrust1_sub <- subset_metabarlist(NC_BerryCrust1, 
                                        table = "pcrs", 
                                        indices = NC_BerryCrust1$pcrs$nb_reads>0 & (
                                          is.na(NC_BerryCrust$pcrs$control_type) |
                                            NC_BerryCrust$pcrs$control_type=="positive"))

# First visualization
comp2 = pcr_within_between(NC_BerryCrust1_sub)
check_pcr_thresh(comp2)

#after removing the non-target and BerryCrust contamination motus, the distance between samples and within samples stayed the same but the density of markers was higher within and lower between.


#now for the "lowering tag-jumps" section:
# Define a vector of thresholds to test
thresholds <- c(0,1e-4,1e-3, 1e-2, 2e-2,3e-2, 5e-2) 

# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(NC_BerryCrust,x, method = "substract"))
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
tmp$controls <- NC_BerryCrust$pcrs$control_type[match(tmp$sample, rownames(NC_BerryCrust$pcrs))]
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
#Removed 34 rows containing non-finite values (`stat_boxplot()`).


#again using the tutorial code which uses cut instead of "substract"

# Define a vector of thresholds to test
thresholds <- c(0,1e-4,1e-3, 1e-2, 2e-2,3e-2, 5e-2) 


# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(NC_BerryCrust,x))
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
tmp$controls <- NC_BerryCrust$pcrs$control_type[match(tmp$sample, rownames(NC_BerryCrust$pcrs))]
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

#Removed 28 rows containing non-finite values (`stat_boxplot()`).

#Did this to the semi-cleaned NC_BerryCrust1 metabarlist to get the MFneglevel of contamination out of the BCneg samples.

# Define a vector of thresholds to test
thresholds <- c(0,1e-4,1e-3, 1e-2, 2e-2,3e-2, 5e-2) 

# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(NC_BerryCrust1,x))
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
tmp$controls <- NC_BerryCrust1$pcrs$control_type[match(tmp$sample, rownames(NC_BerryCrust1$pcrs))]
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

#Removed 582 rows containing non-finite values (`stat_boxplot()`).

#Yay, this looks awesome threshold 0.01 removes all MFneg motus and most of the BCneg motus in the cleaned dataset.

#Summarizing noise in the semi-cleaned dataset

# Create a table of MOTUs quality criteria 
# noise is identified as FALSE in soil_euk, the "!" transforms it to TRUE
motus.qual <- !NC_BerryCrust1$motus[,c("not_a_pcr_conta", "target_taxon", "not_degraded")]
colnames(motus.qual) <- c("pcr_conta", "untargeted_taxon", "degraded_seq")

# Proportion of MOTUs potentially artifactual (TRUE) based on the criteria used
prop.table(table(apply(motus.qual, 1, sum) > 0))

# Corresponding proportion of artifactual reads (TRUE)
prop.table(xtabs(NC_BerryCrust1$motus$count~apply(motus.qual, 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(motus.qual, 2, sum) / nrow(motus.qual)
apply(motus.qual, 2, function(x) sum(NC_BerryCrust1$motus$count[x])/sum(NC_BerryCrust1$motus$count))

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
thresholds <- c(0,1e-4,1e-3, 1e-2, 2e-2,3e-2, 5e-2) 

# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(NC_BerryCrust,x))
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
tmp$controls <- NC_BerryCrust$pcrs$control_type[match(tmp$sample, rownames(NC_BerryCrust$pcrs))]
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

#Removed 28 rows containing non-finite values (`stat_boxplot()`).


# Create a table of MOTUs quality criteria 
# noise is identified as FALSE in NC_BerryCrust, the "!" transforms it to TRUE
motus.qual <- !NC_BerryCrust$motus[,c("not_a_pcr_conta", "target_taxon", "not_degraded")]
colnames(motus.qual) <- c("pcr_conta", "untargeted_taxon", "degraded_seq")

# Proportion of MOTUs potentially artifactual (TRUE) based on the criteria used
prop.table(table(apply(motus.qual, 1, sum) > 0))

# Corresponding proportion of artifactual reads (TRUE)
prop.table(xtabs(NC_BerryCrust$motus$count~apply(motus.qual, 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(motus.qual, 2, sum) / nrow(motus.qual)
apply(motus.qual, 2, function(x) sum(NC_BerryCrust$motus$count[x])/sum(NC_BerryCrust$motus$count))

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

#after cleaning out the non-target taxa and BC neg conta

# Create a table of pcrs quality criteria 
# noise is identified as FALSE in NC_BerryCrust1, the "!" transforms it to TRUE
pcrs.qual <- !NC_BerryCrust1$pcrs[,c("low_contamination_level", "seqdepth_ok")]#, "replicating_pcr")]
colnames(pcrs.qual) <- c("high_contamination_level", "low_seqdepth")#, "outliers")

# Proportion of pcrs potentially artifactual (TRUE) based on the criteria used
# excluding controls
prop.table(table(apply(pcrs.qual[NC_BerryCrust1$pcrs$type=="sample",], 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(pcrs.qual[NC_BerryCrust1$pcrs$type=="sample",], 2, sum) / nrow(pcrs.qual[NC_BerryCrust1$pcrs$type=="sample",])

tmp.pcrs <- 
  apply(sapply(1:ncol(pcrs.qual), function(x) {
    ifelse(pcrs.qual[NC_BerryCrust1$pcrs$type=="sample",x]==T, 
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

#for the full datase

# Create a table of pcrs quality criteria 
# noise is identified as FALSE in NC_BerryCrust, the "!" transforms it to TRUE
pcrs.qual <- !NC_BerryCrust$pcrs[,c("low_contamination_level", "seqdepth_ok")]#, "replicating_pcr")]
colnames(pcrs.qual) <- c("high_contamination_level", "low_seqdepth")#, "outliers")

# Proportion of pcrs potentially artifactual (TRUE) based on the criteria used
# excluding controls
prop.table(table(apply(pcrs.qual[NC_BerryCrust$pcrs$type=="sample",], 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(pcrs.qual[NC_BerryCrust$pcrs$type=="sample",], 2, sum) / nrow(pcrs.qual[NC_BerryCrust$pcrs$type=="sample",])

tmp.pcrs <- 
  apply(sapply(1:ncol(pcrs.qual), function(x) {
    ifelse(pcrs.qual[NC_BerryCrust$pcrs$type=="sample",x]==T, 
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
tmp <- tests[["t_0.01"]]

# Subset on MOTUs: we keep motus that are defined as TRUE following the 
# three criteria below (sum of three TRUE is equal to 3 with the rowSums function)
tmp <- subset_metabarlist(tmp, "motus", 
                          indices = rowSums(tmp$motus[,c("not_a_pcr_conta", "target_taxon")]) == 2)#,
                                                         #"not_degraded")]) == 3)
summary_metabarlist(tmp)

# Subset on pcrs and keep no controls 
NC_BerryCrust_clean <- subset_metabarlist(tmp, "pcrs", 
                                     indices = rowSums(tmp$pcrs[,c("low_contamination_level", "seqdepth_ok")]) == 2 & 
                                       tmp$pcrs$type == "sample")
summary_metabarlist(NC_BerryCrust_clean)

if(sum(colSums(NC_BerryCrust_clean$reads)==0)>0){print("empty motus present")}
if(sum(colSums(NC_BerryCrust_clean$reads)==0)>0){print("empty pcrs present")}

NC_BerryCrust_clean$motus$count = colSums(NC_BerryCrust_clean$reads)
NC_BerryCrust_clean$pcrs$nb_reads_postmetabaR = rowSums(NC_BerryCrust_clean$reads)
NC_BerryCrust_clean$pcrs$nb_motus_postmetabaR = rowSums(ifelse(NC_BerryCrust_clean$reads>0, T, F))

check <- reshape2::melt(NC_BerryCrust_clean$pcrs[,c("nb_reads", "nb_reads_postmetabaR", 
                                     "nb_motus", "nb_motus_postmetabaR")])
check$type <- ifelse(grepl("motus", check$variable), "richness", "abundance")

ggplot(data = check, aes(x = variable, y = value)) +
  geom_boxplot( color = "darkgrey") +
  geom_jitter(alpha=0.1, color = "darkgrey") +
  theme_bw() +
  facet_wrap(~type, scales = "free", ncol = 5) +
  theme(axis.text.x = element_text(angle=45, h=1))


# Compute the number of reads per pcr
NC_BerryCrust_clean$pcrs$nb_reads <- rowSums(NC_BerryCrust_clean$reads)

# Compute the number of motus per pcr
NC_BerryCrust_clean$pcrs$nb_motus <- rowSums(NC_BerryCrust_clean$reads>0)

# Load requested package for plotting
library(ggplot2)
library(reshape2)
library(tidyverse)

# Create an input table (named check2) for ggplot of 3 columns: 
#  (i) control type 
#  (ii) a vector indicated whether it corresponds to nb_reads or nb_motus, 
#  (iii) the corresponding values.

check2 <- reshape2::melt(NC_BerryCrust_clean$pcrs[,c("control_type", "nb_reads", "nb_motus")])

ggplot(data <- check2, aes(x=control_type, y=value, color=control_type)) + 
  geom_boxplot() + theme_bw() + 
  geom_jitter(alpha=0.2) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable, scales = "free_y") + 
  theme(axis.text.x = element_text(angle=45, h=1))

# Using the nb_reads and nb_motus defined previously in the NC_BerryCrust_clean$pcrs table

ggplot(NC_BerryCrust_clean$pcrs, aes(x=nb_reads, y=nb_motus, color = control_type)) + 
  geom_point() + theme_bw() + 
  scale_y_log10() + scale_x_log10() + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey")

ggpcrplate(NC_BerryCrust_clean, FUN = function(m){rowSums(m$reads)}, legend_title = "# of reads per PCR")

#check if this changes measures of beta diversity
# Get row data only for samples
tmp <- subset_metabarlist(NC_BerryCrust, table = "pcrs",
                          indices = NC_BerryCrust$pcrs$type == "sample" &
                            NC_BerryCrust$pcrs$nb_reads>0)

# Add sample biological information for checks
tmp$pcrs$material
#already exists so I don't need to grab a field from Samples for this 
tmp$pcrs$date_collected <- tmp$samples$`Date Collected`[match(tmp$pcrs$sample_id, rownames(tmp$samples))]
tmp$pcrs$date_collected
tmp$pcrs$species <- tmp$samples$`Species (common name)`[match(tmp$pcrs$sample_id, rownames(tmp$samples))]
tmp$pcrs$species

#NC_BerryCrust_clean$pcrs$material #already exists
NC_BerryCrust_clean$pcrs$date_collected <- NC_BerryCrust_clean$samples$`Date Collected`[match(NC_BerryCrust_clean$pcrs$sample_id,
                                                                                      rownames(NC_BerryCrust_clean$samples))]
NC_BerryCrust_clean$pcrs$species <-
  NC_BerryCrust_clean$samples$`Species (common name)`[match(NC_BerryCrust_clean$pcrs$sample_id,
                                                            rownames(NC_BerryCrust_clean$samples))]

#make tmp3 which is NC_BerryCrust_clean without any pcrs with 0 reads.
tmp3 <- subset_metabarlist(NC_BerryCrust_clean, table = "pcrs",
                          indices = #NC_BerryCrust_clean$pcrs$type == "sample" &
                            NC_BerryCrust_clean$pcrs$nb_reads>0)

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
                         #NC_BerryCrust_clean$pcrs$date_collected,
                         sep = " | "))
#Error in vegdist(all, "bray") : 
#missing values are not allowed with argument 'na.rm = FALSE'
#fixed this by adding tmp3 as no empty reads


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

#I tried changing the tag-jump filtering to 0.001 and that added more samples back in, and changing the cutoff back to the pipeline 97% didn't do much.
  
  #in the end I went with 0.01% cutoff for tagjumps and 97% similarity.

#I think it's because the contaminants I blasted were shark sequences (not taxonomically assigned, low confidence etc.) so including the low quality shark motus will make the different species fall out different from each other on the PCA plot.


#For subsetting the data
# Get the samples names from the stomachs plot
stomach_id <- grepl("St", rownames(NC_BerryCrust_clean$pcrs))

# Subset the data
NC_BerryCrust_St <- subset_metabarlist(NC_BerryCrust_clean, table = "pcrs", indices = stomach_id)

# Check results
summary_metabarlist(NC_BerryCrust_St)

#same for fecals
fecal_id<- grepl("NC[1-9]", rownames(NC_BerryCrust_clean$pcrs))
NC_BerryCrust_Fecals <- subset_metabarlist(NC_BerryCrust_clean, table = "pcrs", indices = fecal_id)  
summary_metabarlist(NC_BerryCrust_Fecals)

#now just subset the BerryCrust negative pcr controls
negative_control_id<- grepl("BC", rownames(NC_BerryCrust_clean$pcrs))
NC_BerryCrust_BC_Negatives <- subset_metabarlist(NC_BerryCrust_clean, table = "pcrs", indices = negative_control_id)  
summary_metabarlist(NC_BerryCrust_BC_Negatives)
#hmm 72 motus, 39.5 average!

#now just subset the negative pcr controls
negative_control_id<- grepl("neg", rownames(NC_BerryCrust_clean$pcrs))
NC_BerryCrust_Negatives <- subset_metabarlist(NC_BerryCrust_clean, table = "pcrs", indices = negative_control_id)  
summary_metabarlist(NC_BerryCrust_Negatives)

#View(NC_BerryCrust_Negatives[["motus"]])
#View(NC_BerryCrust_Negatives[["pcrs"]])
#View(NC_BerryCrust_Negatives[["reads"]])
#YAYAYAYAYAY zero motus or reads in the negatives.


#write results to 4 csvs
write.csv(NC_BerryCrust_clean$reads, "NC_BerryCrust_obi3_lulu_metabar_reads.csv")
write.csv(NC_BerryCrust_clean$motus, "NC_BerryCrust_obi3_lulu_metabar_motus.csv")
write.csv(NC_BerryCrust_clean$pcrs, "NC_BerryCrust_obi3_lulu_metabar_pcrs.csv")
write.csv(NC_BerryCrust_clean$samples, "NC_BerryCrust_obi3_lulu_metabar_samples.csv")

#Make a csv that substitutes motus$SCIENTIFIC_NAME for motu_id and puts it as colnames with pcr$rownames and read_numbers from NC_BerryCrust_clean$reads
NC_BerryCrust_clean$pcrs$species <-
  NC_BerryCrust_clean$samples$`Species (scientific name)`[match(NC_BerryCrust_clean$pcrs$sample_id,
rownames(NC_BerryCrust_clean$samples))]

#maybe don't substitute scientific name for motu id, just add it in?

df<-as.data.frame.array(NC_BerryCrust_clean$reads)
rownames(df)
colnames(df)

#maybe don't substitute scientific name for motu id, just add it in????The below command substitutes.
#library(data.table)
#setnames(df, as.character(rownames(NC_BerryCrust_clean$motus)), as.character(NC_BerryCrust_clean$motus$SCIENTIFIC_NAME))


library(tibble)
library(dplyr)
BerryCrust_clean_pcrs<-as.data.frame(NC_BerryCrust_clean$pcrs)
BerryCrust_clean_pcrs<-tibble::rownames_to_column(.data = BerryCrust_clean_pcrs, var="SampleID")
df1<-tibble::rownames_to_column(.data = df,var = "SampleID")
df2<-left_join(df1,BerryCrust_clean_pcrs, by="SampleID")
df3<-df2 %>% relocate(sample_id, .after = SampleID)
df3 <- df3 %>% 
  rename("sample_id" = "SharkID")

df4<-df3 %>% relocate(species, .after = SharkID)
df4<-df4 %>% relocate(date_collected, .after = species)
df4<-df4 %>% relocate(material, .after = date_collected)


#df4<-df4 %>% relocate(nb_reads_postmetabaR, .after = date_collected)
#df4<-df4 %>% relocate(nb_motus_postmetabaR, .after = nb_reads_postmetabaR)

#df5<-df4[,-37:-50]

df5<-column_to_rownames(.data = df4, var = "SampleID")

df6<- t(df5)
#df6<-left_join(df5, NC_BerryCrust_clean$motus)
write.csv(df5, "NC_BerryCrust_clean_species_by_sample_obi3.csv")







#Phyloseq!
#otu_table - Works on any numeric matrix. You must also specify if the species are rows or columns
#sample_data - Works on any data.frame. The rownames must match the sample names in the otu_table if you plan to combine them as a phyloseq-object
#tax_table - Works on any character matrix. The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
#phyloseq - Takes as argument an otu_table and any unordered list of valid phyloseq components: sample_data, tax_table, phylo, or XStringSet. The tip labels of a phylo-object (tree) must match the OTU names of the otu_table, and similarly, the sequence names of an XStringSet object must match the OTU names of the otu_table.
#merge_phyloseq - Can take any number of phyloseq objects and/or phyloseq components, and attempts to combine them into one larger phyloseq object. This is most-useful for adding separately-imported components to an already-created phyloseq object.


#otu_mat <- otu_mat %>%
#  tibble::column_to_rownames("otu")  already imported from lulu
#otu_mat<- as.matrix(otu_mat)
#clean.otu.df<-as.data.frame(t(NC_BerryCrust_clean$reads))
clean.otu.df<-as.data.frame(t(tmp3$reads))
clean.otu.df<-rownames_to_column(clean.otu.df, var="id")

clean.taxa.df<-rownames_to_column(taxa.df.lulu, var="id")
clean.taxa.df<-inner_join(clean.taxa.df,clean.otu.df, by ="id")
#the command above creates a table with taxonomy and number of times each taxa appears in each sample!
clean.taxa.df<-clean.taxa.df[1:8]
clean.taxa.df <- clean.taxa.df %>% 
  tibble::column_to_rownames("id")
samples.df <- samples.df %>% 
  tibble::column_to_rownames("sample") 

clean.taxa.mat <- as.matrix(clean.taxa.df)
clean.otu.df<-column_to_rownames(clean.otu.df, var="id")
clean.otu.mat<-as.matrix(clean.otu.df)

clean.pcrs.df<-as.data.frame((NC_BerryCrust_clean$pcrs))
clean.pcrs.df<-rownames_to_column(clean.pcrs.df, var = "sample")
clean.samples.metabar.df<-rownames_to_column(as.data.frame(NC_BerryCrust_clean$samples),var = "sample_id")


  samples.df.metabar<-full_join(clean.pcrs.df, clean.samples.metabar.df, by= "sample_id")

  samples.df.metabar<-column_to_rownames(samples.df.metabar, var="sample")
  
  
  
 #Starting Phlyoseq!
BerryCrust.ps<- phyloseq(otu_table(clean.otu.mat, taxa_are_rows = TRUE), 
               sample_data(samples.df.metabar), 
               tax_table(clean.taxa.mat))

BerryCrust.ps
sample_names(BerryCrust.ps)
rank_names(BerryCrust.ps)
sample_variables(BerryCrust.ps)

write.csv(clean.otu.df, "ready_for_phyloseq_otus_no_empty_pcrs.csv")
write.csv(samples.df.metabar, "ready_for_phyloseq_samples.csv")
write.csv(clean.taxa.df, "ready_for_phyloseq_taxa.csv")

BerryCrust.ps
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 148 taxa and 102 samples ]
#sample_data() Sample Data:       [ 102 samples by 50 sample variables ]
#tax_table()   Taxonomy Table:    [ 148 taxa by 7 taxonomic ranks ]
#when it says taxa it means motus

plot_richness(BerryCrust.ps, x="material", measures=c("Shannon", "Simpson"), color="material")


# Transform data to proportions as appropriate for Bray-Curtis distances
BerryCrust.ps.prop <- transform_sample_counts(BerryCrust.ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(BerryCrust.ps.prop, method="NMDS", distance="bray")

#Error in if (autotransform && xam > 50) { : 
#missing value where TRUE/FALSE needed

#solved this error by using tmp3(otus from metabar with no empty pcrs (reads>0)) but then got this warning: Warning message:
#In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data

plot_ordination(BerryCrust.ps.prop, ord.nmds.bray, color="species", title="Bray NMDS")
#what does this mean?! changed color to species.

top20 <- names(sort(taxa_sums(BerryCrust.ps), decreasing=TRUE))[0:20]
ps.top20 <- transform_sample_counts(BerryCrust.ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Type", fill="Species") + facet_wrap(~Type, scales="free_x")
#Warning message:
#In psmelt(physeq) : The sample variables: 
#species
#have been renamed to: 
#  sample_species
#to avoid conflicts with taxonomic rank names.

#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(BerryCrust.ps))
standf = function(x, t=total) round(t * (x / sum(x)))
BerryCrust.ps = transform_sample_counts(BerryCrust.ps, standf)

plot_bar(BerryCrust.ps, fill = "genus", facet_grid = "sample_species")
