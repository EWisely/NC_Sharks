#Eldridge Wisely
#2022-12-21
#Eldridge.Wisely@gmail.com
#This script takes the output fastas of Obitools3 for both the crustacean and fish metabarcoding sequences and creates input files for LULU.

#programs and versions
#obitools version 2
#blast version 2.12.0

#change MOTU identifiers to a short unique identifier with the primer name 
cd 04_LULU_Results

#using obitools (not obitools3)
conda activate obitools
#BerryCrust
obiannotate --seq-rank ../03_Obitools3_Results/paired_NC_crustaceans_results.fasta | obiannotate --set-identifier '"BerryCrust_%07d" % seq_rank' > BerryCrust_named.fasta
#sanity check
grep -c 'Callinectes similis' BerryCrust_named.fasta 
#49

#otutable for lulu
obitab -o -d -n 0 BerryCrust_named.fasta >BerryCrust_named.tab
#sanity check
grep -c 'Callinectes similis' BerryCrust_named.tab 
#62

#make a new tab file with only the info Lulu wants.
#just formatting...
#I put the file into excel and added a row to count the field numbers to include.  I didn't save the changes to that file in excel.
cut -f1,7-56 BerryCrust_named.tab > BerryCrust_named_tab_LULU.txt && sed -i -e 's/MERGED_sample://g' BerryCrust_named_tab_LULU.txt


#matchlist for lulu
#First produce a blastdatabase with the OTUs
obiannotate -C BerryCrust_named.fasta >BerryCrust_named_cleared_tags.fasta
#using BLAST
conda activate blast
makeblastdb -in BerryCrust_named_cleared_tags.fasta -parse_seqids -dbtype nucl
# Then blast the OTUs against the database
blastn -db BerryCrust_named_cleared_tags.fasta -outfmt '6 qseqid sseqid pident' -out BerryCrust_named_matchlist.txt -qcov_hsp_perc 80 -perc_identity 84 -query BerryCrust_named_cleared_tags.fasta


#used BerryCrust_named_matchlist.txt and BerryCrust_named_tab_LULU.txt in lulu

# Now, for MiFish

### MiFish name by primer and seq rank, create inputs for lulu
conda activate obitools
obiannotate --seq-rank ../03_Obitools3_Results/paired_NC_fish_results.fasta | obiannotate --set-identifier '"MiFish_%07d" % seq_rank' > MiFish_named.fasta
#checked the input file:
grep -c 'Micropogonias furnieri' ../03_Obitools3_Results/fish_results.fasta
#4
grep -c 'Micropogonias furnieri' MiFish_named.fasta 
#4
#otutable for lulu
obitab -o -d -n 0 MiFish_named.fasta >MiFish_named.tab
grep -c 'Micropogonias furnieri' MiFish_named.tab 
#4

#make a new tab file with only the info Lulu wants.
#I put the file into excel and added a row to count the field numbers to include MiFish_named_tab.xlsx
cut -f1,7-56 MiFish_named.tab > MiFish_named_tab_LULU.txt && sed -i -e 's/MERGED_sample://g' MiFish_named_tab_LULU.txt
#sed -i -e 's/NA/0/g' MiFish_named_tab_LULU.txt


#matchlist for lulu
#First produce a blastdatabase with the OTUs
obiannotate -C MiFish_named.fasta >MiFish_named_cleared_tags.fasta

conda activate blast
makeblastdb -in MiFish_named_cleared_tags.fasta -parse_seqids -dbtype nucl
# Then blast the OTUs against the database
blastn -db MiFish_named_cleared_tags.fasta -outfmt '6 qseqid sseqid pident' -out MiFish_named_matchlist.txt -qcov_hsp_perc 80 -perc_identity 84 -query MiFish_named_cleared_tags.fasta


#used MiFish_named_matchlist.txt and MiFish_named_tab_LULU.txt in lulu
