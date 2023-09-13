#Eldridge Wisely
#2022-12-21
#Eldridge.Wisely@gmail.com
#This script takes the raw data from the UNC sequencing facility with MiFish primers and BerryCrust primers multiplexed together for each sample, and outputs 4 demultiplexed, cleaned, trimmed, annotated, and concatenated files to use as input for Obitools3 pipeline.

#programs and versions
obitools version 2
cutadapt 4.1 with Python 3.9.2

#Optionally, To check if the cleaning works:
#FastQC v0.11.9
#multiqc, version 1.13

#pre-processing steps:

#copy raw data into 01_Data_Raw directory and unzip it

gunzip ../data/00_Data_Raw/*.gz

#remove extra "x" in filenames put in by the sequencing facility.  I did this by hand.
#for example: NC93x-E12_TGGAGTTG-TGTTCCGT_S94_L001_R2_001 becomes NC93-E12_TGGAGTTG-TGTTCCGT_S94_L001_R2_001

#Make a sample list
cd ../data/00_Data_Raw
for files in *R1_001.fastq; do label=$(echo ${files} | cut -d '-' -f '1-2' | cut -d '.' -f 1); echo $label>>../Sample_list.txt; done


#demultiplex MiFish and BerryCrust reads from each sample using Cutadapt and the BerryCrust and MiFish primer sequences.  This makes two files per forward read and two per reverse read.  
cd ..
#pwd should say you're in the data directory
mkdir 01_demultiplexed_MiFish_and_BerryCrust

while read line; do cutadapt -e 0.15 --no-indels -g MiFish=GTCGGTAAAACTCGTGCCAGC -g BerryCrust=GGGACGATAAGACCCTATA -G MiFishUR=CATAGTGGGGTATCTAATCCCAGTTTG -G BerryCrustR=ATTACGCTGTTATCCCTAAAG -o 01_demultiplexed_MiFish_and_BerryCrust/"$line"_R1_{name}.fastq -p 01_demultiplexed_MiFish_and_BerryCrust/"$line"_R2_{name}.fastq 00_Data_Raw/"$line"*R1_001.fastq 00_Data_Raw/"$line"*R2_001.fastq;done<../Sample_list.txt

#Clean and trim demultiplexed reads with cutadapt to get the adapters and stray primers out.
cd 01_demultiplexed_MiFish_and_BerryCrust/
mkdir cutadapt

#Cleaning and trimming the crustacean reads
while read line; do
cutadapt -q 20 --trim-n --minimum-length 10 -a CTTTAGGGATAACAGCGTAATNNN -A TATAGGGTCTTATCGTCCCNNN -G NNNATTACGCTGTTATCCCTAAAG -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CCGAGCCCACGAGAC -a ATCTCGTATGCCGTCTTCTGCTTG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A GACGCTGCCGACGA -A GTGTAGATCTCGGTGGTCGCCGTATCATT -e .1 -o "cutadapt/""$line""_BerryCrust_R1.cleaned.fastq" -p "cutadapt/""$line""_BerryCrust_R2.cleaned.fastq" "$line""_R1_BerryCrust.fastq" "$line""_R2_BerryCrust.fastq"; done < ../Sample_list.txt

mkdir trimmed

while read line; do
trimmomatic PE -trimlog trimmomatic.log "cutadapt/""$line""_BerryCrust_R1.cleaned.fastq" "cutadapt/""$line""_BerryCrust_R2.cleaned.fastq" "trimmed/""$line""_BerryCrust_R1.trimmed.fastq"  "trimmed/""$line""_BerryCrust_R1.discards.fastq" "trimmed/""$line""_BerryCrust_R2.trimmed.fastq" "trimmed/""$line""_BerryCrust_R2.discards.fastq" -validatePairs ILLUMINACLIP:../Nextera_adapters/NexteraPE-PE.fa:2:30:10:2:true SLIDINGWINDOW:4:15;done < ../../Sample_list.txt

#Cleaning and trimming the fish reads
while read line; do
cutadapt -q 20 --trim-n --minimum-length 10 -a CAAACTGGGATTAGATACCCCACTATG -G NNNNCATAGTGGGGTATCTAATCCCAGTTTG -A GCTGGCACGAGTTTTACCGACNNNN -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CCGAGCCCACGAGAC -a ATCTCGTATGCCGTCTTCTGCTTG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A GACGCTGCCGACGA -A GTGTAGATCTCGGTGGTCGCCGTATCATT -e .1 -o "cutadapt/""$line""_MiFish_R1.cleaned.fastq" -p "cutadapt/""$line""_MiFish_R2.cleaned.fastq" "$line""_R1_MiFish.fastq" "$line""_R2_MiFish.fastq"; done < ../Sample_list.txt

while read line; do
trimmomatic PE -trimlog trimmomatic.log "cutadapt/""$line""_MiFish_R1.cleaned.fastq" "cutadapt/""$line""_MiFish_R2.cleaned.fastq" "trimmed/""$line""_MiFish_R1.trimmed.fastq"  "trimmed/""$line""_MiFish_R1.discards.fastq" "trimmed/""$line""_MiFish_R2.trimmed.fastq" "trimmed/""$line""_MiFish_R2.discards.fastq" -validatePairs ILLUMINACLIP:../Nextera_adapters/NexteraPE-PE.fa:2:30:10:2:true SLIDINGWINDOW:4:15;done < ../../Sample_list.txt

#move cleaned data to new folder called 02_Data_Clean
mv trimmed/*.trimmed.fastq ../02_Data_Clean/

#annotate and concatenate the cleaned reads
cd ../02_Data_Clean

#using Obitools (not Obitools3) to annotate with sample and assay (primer) name for each entry in the Fastq files

#First the MiFish reads
while read line; do
obiannotate -S sample:${line} -S assay:MiFish --length ${line}_MiFish_R1.trimmed.fastq > ${line}_MiFish_R1.trimmed.annotated.fastq; done < ../Sample_list.txt

while read line; do
obiannotate -S sample:${line} -S assay:MiFish --length ${line}_MiFish_R2.trimmed.fastq > ${line}_MiFish_R2.trimmed.annotated.fastq; done < ../Sample_list.txt


cat *_MiFish_R2.trimmed.annotated.fastq > NC-MiFish_annotated_R2.fastq
cat *_MiFish_R1.trimmed.annotated.fastq > NC-MiFish_annotated_R1.fastq

#Now the Berry Crustacean reads
while read line; do
obiannotate -S sample:${line} -S assay:BerryCrust --length ${line}_BerryCrust_R1.trimmed.fastq > ${line}_BerryCrust_R1.trimmed.annotated.fastq; done < ../Sample_list.txt

while read line; do
obiannotate -S sample:${line} -S assay:BerryCrust --length ${line}_BerryCrust_R2.trimmed.fastq > ${line}_BerryCrust_R2.trimmed.annotated.fastq; done < ../Sample_list.txt


cat *_BerryCrust_R2.trimmed.annotated.fastq > NC-BerryCrust_annotated_R2.fastq
cat *_BerryCrust_R1.trimmed.annotated.fastq > NC-BerryCrust_annotated_R1.fastq

#Use NC-MiFish_annotated_R2.fastq, NC-MiFish_annotated_R1.fastq, NC-BerryCrust_annotated_R2.fastq, and NC-BerryCrust_annotated_R1.fastq as input for Obitools3 pipeline

