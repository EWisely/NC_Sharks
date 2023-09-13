#back on my laptop
cd Desktop/Savannah/21_Galapagos_diet/
mkdir 06_post-obi3_prep_for_lulu
mv 21Gal_* 06_post-obi3_prep_for_lulu/
conda activate obitools 
 
### BerryCrust name by primer and seq rank, create inputs for lulu
obiannotate --seq-rank 21Gal_BerryCrust_results.fasta | obiannotate --set-identifier '"BerryCrust_%07d" % seq_rank' > BerryCrust_named.fasta
grep -c 'Leptuca deichmanni' BerryCrust_named.fasta 
#3
#otutable for lulu
obitab -o -d -n 0 BerryCrust_named.fasta >BerryCrust_named.tab
grep -c 'Leptuca deichmanni' BerryCrust_named.tab 
#3
make a new tab file with only the info Lulu wants.
#I put the file into excel and added a row to count the field numbers to include
cut -f1,7-180 BerryCrust_named.tab > BerryCrust_named_tab_LULU.txt && sed -i -e 's/MERGED_sample://g' BerryCrust_named_tab_LULU.txt


#matchlist for lulu
#First produce a blastdatabase with the OTUs
obiannotate -C BerryCrust_named.fasta >BerryCrust_named_cleared_tags.fasta
conda activate blast
makeblastdb -in BerryCrust_named_cleared_tags.fasta -parse_seqids -dbtype nucl
# Then blast the OTUs against the database
blastn -db BerryCrust_named_cleared_tags.fasta -outfmt '6 qseqid sseqid pident' -out BerryCrust_named_matchlist.txt -qcov_hsp_perc 80 -perc_identity 84 -query BerryCrust_named_cleared_tags.fasta

### MiFish name by primer and seq rank, create inputs for lulu
obiannotate --seq-rank 21Gal_MiFish_results.fasta | obiannotate --set-identifier '"MiFish_%07d" % seq_rank' > MiFish_named.fasta
grep -c 'Sphoeroides annulatus' MiFish_named.fasta 
#12
#otutable for lulu
obitab -o -d -n 0 MiFish_named.fasta >MiFish_named.tab
grep -c 'Sphoeroides annulatus' MiFish_named.tab 
#12
#make a new tab file with only the info Lulu wants.
#I put the file into excel and added a row to count the field numbers to include
cut -f1,7-180 MiFish_named.tab > MiFish_named_tab_LULU.txt && sed -i -e 's/MERGED_sample://g' MiFish_named_tab_LULU.txt


#matchlist for lulu
#First produce a blastdatabase with the OTUs
obiannotate -C MiFish_named.fasta >MiFish_named_cleared_tags.fasta
conda activate blast
makeblastdb -in MiFish_named_cleared_tags.fasta -parse_seqids -dbtype nucl
# Then blast the OTUs against the database
blastn -db MiFish_named_cleared_tags.fasta -outfmt '6 qseqid sseqid pident' -out MiFish_named_matchlist.txt -qcov_hsp_perc 80 -perc_identity 84 -query MiFish_named_cleared_tags.fasta

#Put it in R to do lulu
