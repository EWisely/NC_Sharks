#programs and versions
obitools version 2
obitools version 3
cutadapt 4.1 with Python 3.9.2
AdapterRemoval ver. 2.3.1
FastQC v0.11.9
multiqc, version 1.13

pre-processing steps:
#demultiplex MiFish and BerryCrust reads using Cutadapt and the BerryCrust and MiFish primers

cutadapt -e 0.15 --no-indels -g MiFish=GTCGGTAAAACTCGTGCCAGC -g BerryCrust=GGGACGATAAGACCCTATA -G MiFishUR=CATAGTGGGGTATCTAATCCCAGTTTG -G BerryCrustR=ATTACGCTGTTATCCCTAAAG -o St96-St97_R1_{name}.fastq -p St96-St97_R2_{name}.fastq St96-St97-J10_CGAAGAAC-TCGAGAGT_S157_L001_R1_001.fastq St96-St97-J10_CGAAGAAC-TCGAGAGT_S157_L001_R2_001.fastq


#BerryCrust:
while read line; do
cutadapt -q 20 --trim-n --minimum-length 10 -a CTTTAGGGATAACAGCGTAATNNN -A TATAGGGTCTTATCGTCCCNNN -G NNNATTACGCTGTTATCCCTAAAG -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CCGAGCCCACGAGAC -a ATCTCGTATGCCGTCTTCTGCTTG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A GACGCTGCCGACGA -A GTGTAGATCTCGGTGGTCGCCGTATCATT -e .1 -o "$line""_BerryCrust_R1.cleaned.fastq" -p "$line""_BerryCrust_R2.cleaned.fastq" "../NC_Sharks_cutadapt/""$line""_R1_BerryCrust.fastq" "../NC_Sharks_cutadapt/""$line""_R2_BerryCrust.fastq"; done < ../Sample_list.txt

#BerryCrust looks pretty good!  I think that perhaps adding the primer-based cleaning at the beginning instead of the end might be helping.
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




#started with the results of cutadapt and trimmomatic:
    MiFish  /Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/clean_w_cutadapt_dec15-22/trimmed/"$line""_MiFish_R1.trimmed.fastq" 
    BerryCrust  /Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/BC_clean_w_cutadapt_dec15-22/trimmed/"$line""_BerryCrust_R1.trimmed.fastq"


MiFish Obitools3 on Panthera

#first annotate and concatenate the cleaned reads
in folder ~/Desktop/Savannah/Cutadapt/cleaned_w_cutadaptdec15-22/trimmed

on my computer with the old obitools installed:
while read line; do
obiannotate -S sample:${line} -S assay:MiFish --length ${line}_MiFish_R1.trimmed.fastq > ${line}_MiFish_R1.trimmed.annotated.fastq; done < ../../Sample_list.txt

while read line; do
obiannotate -S sample:${line} -S assay:MiFish --length ${line}_MiFish_R2.trimmed.fastq > ${line}_MiFish_R2.trimmed.annotated.fastq; done < ../../Sample_list.txt


cat *_MiFish_R2.trimmed.annotated.fastq > NC-MiFish_annotated_R2.fastq
cat *_MiFish_R1.trimmed.annotated.fastq > NC-MiFish_annotated_R1.fastq



#Obitools3 processing of Fish and crustaceans for NC_Sharks project
obi import ../NC-BerryCrust_annotated_R1.fastq crustaceans/reads1
obi import ../NC-BerryCrust_annotated_R2.fastq crustaceans/reads2

obi alignpairedend -R crustaceans/reads2 crustaceans/reads1 crustaceans/aligned_reads

obi stats -a score_norm crustaceans/aligned_reads 
#mean_score_norm	count	total
#         0.897827	6093577	6093577	


obi grep -p "sequence['score_norm'] > 0.8" crustaceans/aligned_reads crustaceans/good_sequences
obi uniq -m sample crustaceans/good_sequences crustaceans/dereplicated_sequences
obi annotate -k COUNT -k MERGED_sample crustaceans/dereplicated_sequences crustaceans/cleaned_metadata_sequences
obi grep -p "len(sequence)>=80 and sequence['COUNT']>=10" crustaceans/cleaned_metadata_sequences crustaceans/denoised_sequences
obi clean -s MERGED_sample -r 0.05 -H crustaceans/denoised_sequences crustaceans/cleaned_sequences_d1_r.5

#time for the reference database


wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

obi import --taxdump taxdump.tar.gz fish/taxonomy/my_tax
obi import --taxdump taxdump.tar.gz crustaceans/taxonomy/my_tax

# import the EMBL invertebrate database and let obi3 make its own db

wget -nH --cut-dirs=6 -A 'STD_INV*.dat.gz' -R 'STD_HUM*.dat.gz','STD_ENV*.dat.gz' -m -np ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/std/
#started at 2:30 PM 12-16-22 on panthera

mkdir EMBL_vrt
cd EMBL_vrt
wget -nH --cut-dirs=6 -A 'STD_VRT*.dat.gz' -R 'STD_HUM*.dat.gz','STD_ENV*.dat.gz' -m -np ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/std/
#started at 8:45AM 21-17-22 on panthera


#BerryCrust
obi import --embl ../EMBL crustaceans/EMBL_inv_refs
obi ecopcr -e 3 -l 100 -L 260 -F GGGACGATAAGACCCTATA -R ATTACGCTGTTATCCCTAAAG --taxonomy crustaceans/taxonomy/my_tax crustaceans/EMBL_inv_refs crustaceans/Berry_ecopcr_refs
obi grep --require-rank=species --require-rank=genus --require-rank=family --taxonomy crustaceans/taxonomy/my_tax crustaceans/Berry_ecopcr_refs crustaceans/Berry_refs_clean

obi build_ref_db -t 0.97 --taxonomy crustaceans/taxonomy/my_tax crustaceans/Berry_refs_clean crustaceans/Berry_refs_clean_97
obi ecotag -m 0.97 --taxonomy crustaceans/taxonomy/my_tax -R crustaceans/Berry_refs_clean_97 crustaceans/cleaned_sequences_d1_r.5 crustaceans/EMBL_assigned_sequences

obi stats -c SCIENTIFIC_NAME crustaceans/EMBL_assigned_sequences



obi align -t 0.95 crustaceans/EMBL_assigned_sequences crustaceans/EMBL_aligned_assigned_sequences

#obi history crustaceans

#obi history -d crustaceans > crustaceans.dot
obi history -d crustaceans/cleaned_sequences_d1_r.5 > crustaceans_one_view.dot

#dot -Tx11 crustaceans.dot

#dot -Tpng crustaceans_one_view.dot -o crustaceans.png
#open crustaceans.png &

obi export --fasta-output crustaceans/EMBL_assigned_sequences -o crustaceans_results.fasta

obi export --tab-output crustaceans/EMBL_aligned_assigned_sequences > crustaceans_results.csv




#MiFish
obi import --embl ../EMBL_vrt fish/EMBL_vrt
obi ecopcr -e 3 -l 100 -L 260 -F GTCGGTAAAACTCGTGCCAGC -R CATAGTGGGGTATCTAATCCCAGTTTG --taxonomy fish/taxonomy/my_tax fish/EMBL_vrt fish/Miya_ecopcr_refs
obi grep --require-rank=species --require-rank=genus --require-rank=family --taxonomy fish/taxonomy/my_tax fish/Miya_ecopcr_refs fish/Miya_refs_clean

obi build_ref_db -t 0.97 --taxonomy fish/taxonomy/my_tax fish/Miya_refs_clean fish/Miya_refs_clean_97


obi import ../NC-MiFish_annotated_R2.fastq fish/reads2
obi import ../NC-MiFish_annotated_R1.fastq fish/reads1

obi alignpairedend -R fish/reads2 fish/reads1 fish/aligned_reads

#mean_score_norm	count	total
#          0.931502	3942150	3942150	


obi grep -p "sequence['score_norm'] > 0.8" fish/aligned_reads fish/good_sequences
obi uniq -m sample fish/good_sequences fish/dereplicated_sequences
obi annotate -k COUNT -k MERGED_sample fish/dereplicated_sequences fish/cleaned_metadata_sequences
obi grep -p "len(sequence)>=80 and sequence['COUNT']>=10" fish/cleaned_metadata_sequences fish/denoised_sequences
obi clean -s MERGED_sample -r 0.05 -H fish/denoised_sequences fish/cleaned_sequences_d1_r.5


obi ecotag -m 0.97 --taxonomy fish/taxonomy/my_tax -R fish/Miya_refs_clean_97 fish/cleaned_sequences_d1_r.5 fish/EMBL_assigned_sequences

obi stats -c SCIENTIFIC_NAME fish/EMBL_assigned_sequences

obi align -t 0.95 fish/EMBL_assigned_sequences fish/EMBL_aligned_assigned_sequences
obi history -d fish/cleaned_sequences_d1_r.5 > fish_one_view.dot

obi export --fasta-output fish/EMBL_assigned_sequences -o fish_results.fasta

obi export --tab-output fish/EMBL_aligned_assigned_sequences > fish_results.csv



#all background dead-ends and failed attempts and roads-not-taken recorded in ~/Desktop/Savannah/NC_Sharks_data/NC Sharks Data Analysis Steps and Commands.txt