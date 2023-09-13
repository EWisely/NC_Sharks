#Obitools3 processing of Fish and crustaceans for NC_Sharks project

cd 03_Obitools3_Results
#activate obitools3 environment, then:
#Following the tutorial on https://git.metabarcoding.org/obitools/obitools3/-/wikis/Wolf-tutorial-with-the-OBITools3

#Import reads
obi import ../02_Data_Clean/NC-BerryCrust_annotated_R1.fastq crustaceans/reads1
obi import ../02_Data_Clean/NC-BerryCrust_annotated_R2.fastq crustaceans/reads2

#Align forward and reverse reads into full sequences
obi alignpairedend -R crustaceans/reads2 crustaceans/reads1 crustaceans/aligned_reads

#Inspect the alignment scores from the last step
obi stats -a score_norm crustaceans/aligned_reads 
#mean_score_norm	count	total
#         0.897827	6093577	6093577	

#Get just the sequences with over 80% alignment score
obi grep -p "sequence['score_norm'] > 0.8" crustaceans/aligned_reads crustaceans/good_sequences
#Dereplicate the good sequences and count how many replicates of each sequence were in each sample
obi uniq -m sample crustaceans/good_sequences crustaceans/dereplicated_sequences
#clean up the annotations in the dereplicated sequences, keeping only two fields from the previous steps: COUNT, and MERGED_sample
obi annotate -k COUNT -k MERGED_sample crustaceans/dereplicated_sequences crustaceans/cleaned_metadata_sequences
#Get just the sequences longer than 80 basepairs, and which occur in the entire dataset more than 10 times
obi grep -p "len(sequence)>=80 and sequence['COUNT']>=10" crustaceans/cleaned_metadata_sequences crustaceans/denoised_sequences
#Obi clean checks for PCR artifacts and groups sequences based on sequence similarity and relative abundance thresholds.  Sequences expected to be real are called "Heads"
obi clean -s MERGED_sample -r 0.05 -H crustaceans/denoised_sequences crustaceans/cleaned_sequences_d1_r.5

#time for the reference database
#get the NCBI taxdump.  This is not a huge file and doesn't take long.
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
obi import --taxdump taxdump.tar.gz fish/taxonomy/my_tax
obi import --taxdump taxdump.tar.gz crustaceans/taxonomy/my_tax

# import the EMBL invertebrate database and let obi3 make its own db
#This does take a long time and lots of disk space.  The command downloads the Invertebrate database from EMBL, and specifically excludes anything from the Human or Environmental Datasets
mkdir EMBL_inv
cd EMBL_inv
wget -nH --cut-dirs=6 -A 'STD_INV*.dat.gz' -R 'STD_HUM*.dat.gz','STD_ENV*.dat.gz' -m -np ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/std/
#started at 2:30 PM 12-16-22



#For the crustaceans 
obi import --embl ../EMBL crustaceans/EMBL_inv_refs
#do ecopcr with Berry's Crustacean primers to pull out all invertebrate sequences that could amplify with the primers we used (allowing 3 errors in the primer binding regions)
obi ecopcr -e 3 -l 100 -L 260 -F GGGACGATAAGACCCTATA -R ATTACGCTGTTATCCCTAAAG --taxonomy crustaceans/taxonomy/my_tax crustaceans/EMBL_inv_refs crustaceans/Berry_ecopcr_refs
#Get only the sequences that have species, genus, and family information associated with them.
obi grep --require-rank=species --require-rank=genus --require-rank=family --taxonomy crustaceans/taxonomy/my_tax crustaceans/Berry_ecopcr_refs crustaceans/Berry_refs_clean
#Build ecopcr reference database to use in the ecotag taxonomic assignment step, threshold 97% sequence similarity
obi build_ref_db -t 0.97 --taxonomy crustaceans/taxonomy/my_tax crustaceans/Berry_refs_clean crustaceans/Berry_refs_clean_97
#Ecotag all of the BerryCrust assay sequences with their taxonomic information if they're at least 97% similar to something in the database.  Runs Lowest Common Ancestor (LCA) analysis to determine if sequences can only be ID'ed to genus or family level.
obi ecotag -m 0.97 --taxonomy crustaceans/taxonomy/my_tax -R crustaceans/Berry_refs_clean_97 crustaceans/cleaned_sequences_d1_r.5 crustaceans/EMBL_assigned_sequences
#View results!
obi stats -c SCIENTIFIC_NAME crustaceans/EMBL_assigned_sequences
#results show up here for the full dataset

#Optionally align the sequences and look at the history files
obi align -t 0.95 crustaceans/EMBL_assigned_sequences crustaceans/EMBL_aligned_assigned_sequences

#obi history crustaceans

#obi history -d crustaceans > crustaceans.dot
obi history -d crustaceans/cleaned_sequences_d1_r.5 > crustaceans_one_view.dot

#dot -Tx11 crustaceans.dot

#dot -Tpng crustaceans_one_view.dot -o crustaceans.png
#open crustaceans.png &

#Export crustaceans_results.fasta
obi export --fasta-output crustaceans/EMBL_assigned_sequences -o crustaceans_results.fasta

#obi export --tab-output crustaceans/EMBL_aligned_assigned_sequences > crustaceans_results.csv




#MiFish

#Get the EMBL vertebrate database, not including human sequences or environmental sequences.
mkdir EMBL_vrt
cd EMBL_vrt
wget -nH --cut-dirs=6 -A 'STD_VRT*.dat.gz' -R 'STD_HUM*.dat.gz','STD_ENV*.dat.gz' -m -np ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/std/
#started at 8:45AM 12-17-22

#go back to 03_Obitools3_Results directory
cd ..
#import database into the obitools3 DMS called fish
obi import --embl EMBL_vrt fish/EMBL_vrt
#in-silico PCR with MiFish primers allowing up to 3 mismatches.  
#This produces a file with all of the vertebrates in the EMBL database that could amplify with the MiFish primers.
obi ecopcr -e 3 -l 100 -L 260 -F GTCGGTAAAACTCGTGCCAGC -R CATAGTGGGGTATCTAATCCCAGTTTG --taxonomy fish/taxonomy/my_tax fish/EMBL_vrt fish/Miya_ecopcr_refs
#only use the sequences that have species, genus, and family information attached.
obi grep --require-rank=species --require-rank=genus --require-rank=family --taxonomy fish/taxonomy/my_tax fish/Miya_ecopcr_refs fish/Miya_refs_clean
#build the ecotag reference database, threshold 97%
obi build_ref_db -t 0.97 --taxonomy fish/taxonomy/my_tax fish/Miya_refs_clean fish/Miya_refs_clean_97

#Import the MiFish forward and reverse reads into the fish DMS
obi import ../02_Data_Clean/NC-MiFish_annotated_R2.fastq fish/reads2
obi import ../02_Data_Clean/NC-MiFish_annotated_R1.fastq fish/reads1

#align forward and reverse reads
obi alignpairedend -R fish/reads2 fish/reads1 fish/aligned_reads

#mean_score_norm	count	total
#          0.931502	3942150	3942150	

#use only the ones with better than 80% alignment scores
obi grep -p "sequence['score_norm'] > 0.8" fish/aligned_reads fish/good_sequences
#dereplicate and count
obi uniq -m sample fish/good_sequences fish/dereplicated_sequences
#clean up the annotations of the dereplicated sequences
obi annotate -k COUNT -k MERGED_sample fish/dereplicated_sequences fish/cleaned_metadata_sequences
#use only the ones longer than 80 basepairs and with more than 10 occurences in the dataset
obi grep -p "len(sequence)>=80 and sequence['COUNT']>=10" fish/cleaned_metadata_sequences fish/denoised_sequences
#obi clean to filter out suspected PCR artifacts and get the final list of MiFish MOTUS
obi clean -s MERGED_sample -r 0.05 -H fish/denoised_sequences fish/cleaned_sequences_d1_r.5

#assign taxonomy to those that match a sequence in the database with more than 97% sequence similarity and pass LCA filtering
obi ecotag -m 0.97 --taxonomy fish/taxonomy/my_tax -R fish/Miya_refs_clean_97 fish/cleaned_sequences_d1_r.5 fish/EMBL_assigned_sequences
#look at the results!
obi stats -c SCIENTIFIC_NAME fish/EMBL_assigned_sequences

obi align -t 0.95 fish/EMBL_assigned_sequences fish/EMBL_aligned_assigned_sequences
obi history -d fish/cleaned_sequences_d1_r.5 > fish_one_view.dot

#export the results out of the obitools3 DMS
obi export --fasta-output fish/EMBL_assigned_sequences -o fish_results.fasta

#obi export --tab-output fish/EMBL_aligned_assigned_sequences > fish_results.csv
