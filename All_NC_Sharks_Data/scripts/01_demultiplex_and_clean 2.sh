#Eldridge Wisely
#2022-12-21
#Eldridge.Wisely@gmail.com
#This script takes the raw data from the UNC sequencing facility with MiFish primers and BerryCrust primers multiplexed together for each sample, and outputs 4 demultiplexed, cleaned, trimmed, annotated, and concatenated files to use as input for Obitools3 pipeline.

#programs and versions
#obitools version 2
#cutadapt 4.1 with Python 3.9.2

#Optionally, To check if the cleaning works:
#FastQC v0.11.9
#multiqc, version 1.13

#pre-processing steps:

#copy raw data into 00_Data_Raw directory and untar it

tar -xzvf 00_Paired_Data_Raw.tar.gz

#remove extra "x" in filenames put in by the sequencing facility.  
#for example: NC93x-E12_TGGAGTTG-TGTTCCGT_S94_L001_R2_001 becomes NC93-E12_TGGAGTTG-TGTTCCGT_S94_L001_R2_001

for file_name in *.gz
do 
  new_file_name=$(sed 's/[x]*\-/-/g' <<< "$file_name");
  mv "$file_name" "$new_file_name";
done

#Make a sample list
cd ../data/00_Data_Raw
for files in *R1_001.fastq.gz; do label=$(echo ${files} | cut -d '_' -f '1-2' | cut -d '.' -f 1); echo $label>>../Sample_list.txt; done


#demultiplex MiFish and BerryCrust reads from each sample using Cutadapt and the BerryCrust and MiFish primer sequences.  This makes two files per forward read and two per reverse read.  
cd ..
#pwd should say you're in the data directory
mkdir 01_demultiplexed_MiFish_and_BerryCrust

conda activate qc

while read line; do cutadapt -e 0.15 --no-indels -g MiFish=GTCGGTAAAACTCGTGCCAGC -g BerryCrust=GGGACGATAAGACCCTATA -G MiFishUR=CATAGTGGGGTATCTAATCCCAGTTTG -G BerryCrustR=ATTACGCTGTTATCCCTAAAG -o 01_demultiplexed_MiFish_and_BerryCrust/"$line"_R1_{name}.fastq.gz -p 01_demultiplexed_MiFish_and_BerryCrust/"$line"_R2_{name}.fastq.gz 00_Data_Raw/"$line"*R1_001.fastq.gz 00_Data_Raw/"$line"*R2_001.fastq.gz;done<Sample_list.txt

#Clean and trim demultiplexed reads with cutadapt to get the adapters and stray primers out.
cd 01_demultiplexed_MiFish_and_BerryCrust/
mkdir cutadapt

#Cleaning and trimming the crustacean reads
while read line; do
cutadapt -q 20 --trim-n --minimum-length 10 -a CTTTAGGGATAACAGCGTAATNNN -A TATAGGGTCTTATCGTCCCNNN -G NNNATTACGCTGTTATCCCTAAAG -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CCGAGCCCACGAGAC -a ATCTCGTATGCCGTCTTCTGCTTG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A GACGCTGCCGACGA -A GTGTAGATCTCGGTGGTCGCCGTATCATT -e .1 -o "cutadapt/""$line""_BerryCrust_R1.cleaned.fastq.gz" -p "cutadapt/""$line""_BerryCrust_R2.cleaned.fastq.gz" "$line""_R1_BerryCrust.fastq.gz" "$line""_R2_BerryCrust.fastq.gz"; done < ../Sample_list.txt
#Download the NexteraPE-PE.fa from Trimmomatic's website and put it here: ../Nextera_adapters/NexteraPE-PE.fa  (Also available in this github: data/Nextera_adapters/NexteraPE-PE.fa)

mkdir trimmed

while read line; do
trimmomatic PE -trimlog trimmomatic.log "cutadapt/""$line""_BerryCrust_R1.cleaned.fastq.gz" "cutadapt/""$line""_BerryCrust_R2.cleaned.fastq.gz" "trimmed/""$line""_BerryCrust_R1.trimmed.fastq.gz"  "trimmed/""$line""_BerryCrust_R1.discards.fastq.gz" "trimmed/""$line""_BerryCrust_R2.trimmed.fastq.gz" "trimmed/""$line""_BerryCrust_R2.discards.fastq.gz" -validatePairs ILLUMINACLIP:../Nextera_adapters/NexteraPE-PE.fa:2:30:10:2:true SLIDINGWINDOW:4:15;done < ../Sample_list.txt

#Cleaning and trimming the fish reads
while read line; do
cutadapt -q 20 --trim-n --minimum-length 10 -a CAAACTGGGATTAGATACCCCACTATG -G NNNNCATAGTGGGGTATCTAATCCCAGTTTG -A GCTGGCACGAGTTTTACCGACNNNN -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CCGAGCCCACGAGAC -a ATCTCGTATGCCGTCTTCTGCTTG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A GACGCTGCCGACGA -A GTGTAGATCTCGGTGGTCGCCGTATCATT -e .1 -o "cutadapt/""$line""_MiFish_R1.cleaned.fastq.gz" -p "cutadapt/""$line""_MiFish_R2.cleaned.fastq.gz" "$line""_R1_MiFish.fastq.gz" "$line""_R2_MiFish.fastq.gz"; done < ../Sample_list.txt

while read line; do
trimmomatic PE -trimlog trimmomatic.log "cutadapt/""$line""_MiFish_R1.cleaned.fastq.gz" "cutadapt/""$line""_MiFish_R2.cleaned.fastq.gz" "trimmed/""$line""_MiFish_R1.trimmed.fastq.gz"  "trimmed/""$line""_MiFish_R1.discards.fastq.gz" "trimmed/""$line""_MiFish_R2.trimmed.fastq.gz" "trimmed/""$line""_MiFish_R2.discards.fastq.gz" -validatePairs ILLUMINACLIP:../Nextera_adapters/NexteraPE-PE.fa:2:30:10:2:true SLIDINGWINDOW:4:15;done < ../Sample_list.txt

#move cleaned data to new folder called 02_Data_Clean
mkdir ../02_Data_Clean
cp trimmed/*.trimmed.fastq.gz ../02_Data_Clean/

#annotate and concatenate the cleaned reads
cd ../02_Data_Clean

for file_name in *.gz
do 
  new_file_name=$(sed 's/-[^S]*\_B/_B/g' <<< "$file_name");
  mv "$file_name" "$new_file_name";
done

for file_name in *.gz
do 
  new_file_name=$(sed 's/-[^S]*\_M/_M/g' <<< "$file_name");
  mv "$file_name" "$new_file_name";
done

#make new Sample_list.txt without the sequencing tag names 

gunzip *.fastq.gz
conda activate obitools

#using Obitools (not Obitools3) to annotate with sample and assay (primer) name for each entry in the Fastq files


# This next section is for just the paired samples, use the full sample list and the full 00_Raw_Data from ORIGINAL folder for the bigger project.

#First the MiFish reads
while read line; do
obiannotate -S sample:${line} -S assay:MiFish --length ${line}_MiFish_R1.trimmed.fastq > ${line}_MiFish_R1.trimmed.annotated.fastq; done < ../Paired_sample_list.txt

while read line; do
obiannotate -S sample:${line} -S assay:MiFish --length ${line}_MiFish_R2.trimmed.fastq > ${line}_MiFish_R2.trimmed.annotated.fastq; done < ../Paired_sample_list.txt


cat *_MiFish_R2.trimmed.annotated.fastq > Paired_NC-MiFish_annotated_R2.fastq
cat *_MiFish_R1.trimmed.annotated.fastq > Paired_NC-MiFish_annotated_R1.fastq

#Now the Berry Crustacean reads
while read line; do
obiannotate -S sample:${line} -S assay:BerryCrust --length ${line}_BerryCrust_R1.trimmed.fastq > ${line}_BerryCrust_R1.trimmed.annotated.fastq; done < ../Paired_sample_list.txt

while read line; do
obiannotate -S sample:${line} -S assay:BerryCrust --length ${line}_BerryCrust_R2.trimmed.fastq > ${line}_BerryCrust_R2.trimmed.annotated.fastq; done < ../Paired_sample_list.txt


cat *_BerryCrust_R2.trimmed.annotated.fastq > Paired_NC-BerryCrust_annotated_R2.fastq
cat *_BerryCrust_R1.trimmed.annotated.fastq > Paired_NC-BerryCrust_annotated_R1.fastq

#Use NC-MiFish_annotated_R2.fastq, NC-MiFish_annotated_R1.fastq, NC-BerryCrust_annotated_R2.fastq, and NC-BerryCrust_annotated_R1.fastq as input for Obitools3 pipeline for the full dataset
#Use Paired_NC-MiFish_annotated_R2.fastq, Paired_NC-MiFish_annotated_R1.fastq, Paired_NC-BerryCrust_annotated_R2.fastq, and Paired_NC-BerryCrust_annotated_R1.fastq as input for Obitools3 pipeline for the full dataset


