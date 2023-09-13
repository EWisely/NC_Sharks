#programs and versions
obitools version 2
obitools version 3
cutadapt 4.1 with Python 3.9.2
AdapterRemoval ver. 2.3.1
FastQC v0.11.9
multiqc, version 1.13

#pre-processing steps:
#demultiplex MiFish and BerryCrust reads using Cutadapt and the BerryCrust and MiFish primers

#example with one sample:
cutadapt -e 0.15 --no-indels -g MiFish=GTCGGTAAAACTCGTGCCAGC -g BerryCrust=GGGACGATAAGACCCTATA -G MiFishUR=CATAGTGGGGTATCTAATCCCAGTTTG -G BerryCrustR=ATTACGCTGTTATCCCTAAAG -o St96-St97_R1_{name}.fastq -p St96-St97_R2_{name}.fastq St96-St97-J10_CGAAGAAC-TCGAGAGT_S157_L001_R1_001.fastq St96-St97-J10_CGAAGAAC-TCGAGAGT_S157_L001_R2_001.fastq


#BerryCrust with primer-based cleaning
while read line; do
cutadapt -q 20 --trim-n --minimum-length 10 -a CTTTAGGGATAACAGCGTAATNNN -A TATAGGGTCTTATCGTCCCNNN -G NNNATTACGCTGTTATCCCTAAAG -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CCGAGCCCACGAGAC -a ATCTCGTATGCCGTCTTCTGCTTG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A GACGCTGCCGACGA -A GTGTAGATCTCGGTGGTCGCCGTATCATT -e .1 -o "$line""_BerryCrust_R1.cleaned.fastq" -p "$line""_BerryCrust_R2.cleaned.fastq" "../NC_Sharks_cutadapt/""$line""_R1_BerryCrust.fastq" "../NC_Sharks_cutadapt/""$line""_R2_BerryCrust.fastq"; done < ../Sample_list.txt

mkdir trimmed
cd trimmed
while read line; do
trimmomatic PE -trimlog trimmomatic.log "../""$line""_BerryCrust_R1.cleaned.fastq" "../""$line""_BerryCrust_R2.cleaned.fastq" "$line""_BerryCrust_R1.trimmed.fastq"  "$line""_BerryCrust_R1.discards.fastq" "$line""_BerryCrust_R2.trimmed.fastq" "$line""_BerryCrust_R2.discards.fastq" -validatePairs ILLUMINACLIP:../../NexteraPE-PE.fa:2:30:10:2:true SLIDINGWINDOW:4:15;done < ../../Sample_list.txt

#MiFish:
#primer-based cleaning added to the front of the command
while read line; do
cutadapt -q 20 --trim-n --minimum-length 10 -a CAAACTGGGATTAGATACCCCACTATG -G NNNNCATAGTGGGGTATCTAATCCCAGTTTG -A GCTGGCACGAGTTTTACCGACNNNN -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CCGAGCCCACGAGAC -a ATCTCGTATGCCGTCTTCTGCTTG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A GACGCTGCCGACGA -A GTGTAGATCTCGGTGGTCGCCGTATCATT -e .1 -o "$line""_MiFish_R1.cleaned.fastq" -p "$line""_MiFish_R2.cleaned.fastq" "../NC_Sharks_cutadapt/""$line""_R1_MiFish.fastq" "../NC_Sharks_cutadapt/""$line""_R2_MiFish.fastq"; done < ../Sample_list.txt

while read line; do
trimmomatic PE -trimlog trimmomatic.log "../""$line""_MiFish_R1.cleaned.fastq" "../""$line""_MiFish_R2.cleaned.fastq" "$line""_MiFish_R1.trimmed.fastq"  "$line""_MiFish_R1.discards.fastq" "$line""_MiFish_R2.trimmed.fastq" "$line""_MiFish_R2.discards.fastq" -validatePairs ILLUMINACLIP:../../NexteraPE-PE.fa:2:30:10:2:true SLIDINGWINDOW:4:15;done < ../../Sample_list.txt

