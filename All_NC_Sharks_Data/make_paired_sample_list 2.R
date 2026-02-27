library(tidyverse)
all_samples_metadata<-read_delim("data/05_MetabaR_inputs/MetabaR_formatted_samples_NCSharks.txt")
paired_column_added<-read.csv("data/05_MetabaR_inputs/newfilewithpairedcolumn.csv")
paired_samples<-paired_column_added%>%
  filter(Paired=="Yes")
paired_shark_list<-paired_sample_metadata%>%
  select(sample_id)

paired_samples_metadata<-inner_join(all_samples_metadata,paired_shark_list)
view(paired_samples_metadata)

write_delim(paired_samples_metadata, file="../Paired_only/Paired_MetabaR_formatted_samples_NCSharks.txt")

paired_PCR_list<-paired_samples%>%
  select(X)%>%
  rename(sample=X)


all_pcrs_metadata<-read_delim("data/05_MetabaR_inputs/MiFish_Metabar_formatted_pcrs_pureobi3_lulu_noSt9.txt")
control_PCRs_metadata<-all_pcrs_metadata%>%
  filter(type=="control")

paired_pcrs_metadata<-inner_join(all_pcrs_metadata,paired_PCR_list)
paired_pcrs_with_controls<-full_join(paired_pcrs_metadata,control_PCRs_metadata)

paired_pcrs_with_controls_list<-paired_pcrs_with_controls%>%
  select(sample)

write_delim(paired_pcrs_with_controls, file="../Paired_only/Paired_Metabar_formatted_PCRs.txt")

write_delim(paired_pcrs_with_controls_list, file = "../Paired_only/Paired_sample_list.txt",col_names = F)

#use the following code to copy the original data from the original data file to the paired only directory.

#while read line; do cp -- $line[x-]*.fastq.gz ../following_obitools3_method/Paired_only/00_Data_Raw/;done<../following_obitools3_method/Paired_only/Paired_sample_list.txt

#upload this to Figshare

