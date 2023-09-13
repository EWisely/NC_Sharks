eDNA/fDNA commands (used for NC Sharks Oct 2022)

starting from demultiplexed (from the sequencing facility) reads, I will try to treat them as un-demultiplexed because I want to use 
eDNA flow for quality control as well and submitting the demultiplexed reads presents some questions:
 (do I need to keep the primer sequences or not? How do I put the fastq files in usearch format?)
 
 OK, I will also try pulling the demuxed with cutadapt reads directly from google drive:
 https://drive.google.com/drive/folders/107JqMeXUgcLfnTYXuM8sNkqrIrXpmJiQ?usp=sharing

command from : https://medium.com/@acpanjan/download-google-drive-files-using-wget-3c2c025a8b99

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=107JqMeXUgcLfnTYXuM8sNkqrIrXpmJiQ' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=107JqMeXUgcLfnTYXuM8sNkqrIrXpmJiQ" -O FILENAME && rm -rf /tmp/cookies.txt
started to work but I think it failed because I couldn't make the file public access enough, just anyone with link.  Because it's U of A google drive.

oh well, it's only 2.5GB, I can download it to my computer and then copy it to leopardus

my computer is full.  ugh.

PARAMETERS for default eDNAFlow run:
from eDNAFlow nextflow.config file on github:
params {
reads   = "*_{R1,R2}.fastq"
barcode = "*_bc.txt"
minQuality = "20"
minAlignLeng = "12"
minLen = "50"
primer_mismatch = "2"
minsize = "8"
maxTarSeq = "10"
perc_identity = "95"
evalue = "1e-3"
qcov = "100"
lulu = "lulu.R"
mode = "usearch32"  
usearch64 = ""   
blast_db = ""
custom_db = ""
blast_task = "blastn"
publish
_dir_mode = "symlink" 
bindDir = ""
singularityDir = ""

lca_script = "LCA_taxonomyAssignment_scripts/runAssign_collapsedTaxonomy.py"
zotuTable = ""
blastFile = ""
lca_qcov = "100"
lca_pid = "97"
lca_diff = "1"
lca_output = "lca_taxAssigned_results"

following the eDNAFlow pipeline by hand until I can drop it in after the demultiplexing step:
ran fastqc and multiqc on the demultiplexed files that Savannah sent back (leopardus qc environment)
next: adapterremoval
AdapterRemoval --threads ${task.cpus} --file1 ${read[0]} --file2 ${read[1]} \
                    --collapse --trimns --trimqualities \
                    --minquality $minQuality \
                    --minalignmentlength $minAlignLeng \
                    --basename ${sample_id}
    mv ${sample_id}.collapsed ${sample_id}_QF.fastq  

while read line; do    
AdapterRemoval --threads 10 --file1 "../MiFish/""$line""_R1_MiFish.fastq" --file2 "../MiFish/""$line""_R2_MiFish.fastq" \
                    --collapse --trimns --trimqualities \
                    --minquality 20 \
                    --minalignmentlength 12 \
                    --basename "$line""_MiFish"; done < ../Sample_list.txt
                    
while read line; do    
AdapterRemoval --threads 10 --file1 "../BerryCrust/""$line""_R1_BerryCrust.fastq" --file2 "../BerryCrust/""$line""_R2_BerryCrust.fastq" \
                    --collapse --trimns --trimqualities \
                    --minquality 20 \
                    --minalignmentlength 12 \
                    --basename "$line""_BerryCrust"; done < ../Sample_list.txt
                    
Worked!  After converting the Sample_list.txt using the dos2unix command on leopardus

mkdir ../01b_fastqc_merged_filtered
fastqc *.collapsed -o ../01b_fastqc_merged_filtered/
cd ../01b_fastqc_merged_filtered/
multiqc .

I think they all need adapter removal done because the multiqc still says many of them have Nextera transposase adapter contamination.  The lengths look really good though.
I checked some of the files of adapter contamination concern from fastqc/multiqc and both the MiFish and BerryCrust reads still have the forward PCR1 primer adapter TCGTCGGCAGCGTC and often followed by the 3' end of the illumina adapter AGATGTGTATAAGAGACAG hanging out on the ends of the reads.

I'm going to try cutadapt -b TCGTCGGCAGCGTC -o output.fastq inputfile.collapsed
while read line; do    
cutadapt -b TCGTCGGCAGCGTC -o "../01c_cutadapt_merged_qualityfiltered/""$line""_BerryCrust.collapsed.AR" "$line""_BerryCrust.collapsed" ; done < ../Sample_list.txt
worked!

Now the same thing for MiFish
while read line; do    
cutadapt -b TCGTCGGCAGCGTC -o "../01c_cutadapt_merged_qualityfiltered/""$line""_MiFish.collapsed.AR" "$line""_MiFish.collapsed" ; done < ../Sample_list.txt
 

Now I'll fastqc and multiqc these results to make sure they're clean
cd ../01c_cutadapt_merged_qualityfiltered/
mkdir ../01d_fastqc_merged_filtered_adapter_removed
fastqc * -o ../01d_fastqc_merged_filtered_adapter_removed/
cd ../01d_fastqc_merged_filtered_adapter_removed
multiqc .

looks like there's still adapter in a lot of the BerryCrust ones and there are some BerryCrust empty output folders from when I accidentally forgot the _before BerryCrust.  

ok, I got rid of the empty files because of the typo.
Now, to do the above commands with the second half of the PCR1 forward primer AGATGTGTATAAGAGACAG

while read line; do    
cutadapt -b AGATGTGTATAAGAGACAG -o "$line""_BerryCrust.collapsed.AR1" "first_adapter_TCGTCGGCAGCGTC/""$line""_BerryCrust.collapsed.AR" ; done < ../Sample_list.txt

while read line; do    
cutadapt -b AGATGTGTATAAGAGACAG -o "$line""_MiFish.collapsed.AR1" "first_adapter_TCGTCGGCAGCGTC/""$line""_MiFish.collapsed.AR" ; done < ../Sample_list.txt

#in 01c_cutadapt_merged_qualityfiltered
mkdir ../01e_fastqc_second_adapter_removed_AR1
fastqc *.AR1 -o ../01e_fastqc_second_adapter_removed_AR1
cd ../01e_fastqc_second_adapter_removed_AR1
multiqc .

nope!  found that the reverse complement of the same adapters from PCR1 were still in lots of sequences CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
and the same with the reverse complement of the adapter used in the second PCR ATCTCGTATGCCGTCTTCTGCTTG
and the adapter used in the forward primer of PCR2 AATGATACGGCGACCACCGAGATCTACAC

to remove them all at once from the original.collapsed files:
deleted the previous attempts to save disk space
rm -r 01e_fastqc_second_adapter_removed_AR1
rm -r 01d_fastqc_merged_filtered_adapter_removed
rm 01c_cutadapt_merged_qualityfiltered/*

cd 01_merged_qualityfiltered/
while read line; do    
cutadapt --revcomp --times 10 -b TCGTCGGCAGCGTC -b AGATGTGTATAAGAGACAG -b GTCTCGTGGGCTCGG -b AATGATACGGCGACCACCGAGATCTACAC -b CAAGCAGAAGACGGCATACGAGAT -o "../01c_cutadapt_merged_qualityfiltered/""$line""_BerryCrust.collapsed.AR" "$line""_BerryCrust.collapsed" ; done < ../Sample_list.txt
worked!

while read line; do    
cutadapt --revcomp --times 10 -b TCGTCGGCAGCGTC -b AGATGTGTATAAGAGACAG -b GTCTCGTGGGCTCGG -b AATGATACGGCGACCACCGAGATCTACAC -b CAAGCAGAAGACGGCATACGAGAT -o "../01c_cutadapt_merged_qualityfiltered/""$line""_MiFish.collapsed.AR" "$line""_MiFish.collapsed" ; done < ../Sample_list.txt
worked



cd ../01c_cutadapt_merged_qualityfiltered/
mkdir ../01d_fastqc_merged_filtered_adapter_removed
fastqc * -o ../01d_fastqc_merged_filtered_adapter_removed/
cd ../01d_fastqc_merged_filtered_adapter_removed
multiqc .

looks like there are still some revcomp adapters especially in the BerryCrust reads
while read line; do    
cutadapt --times 10 -b TCGTCGGCAGCGTC -b AGATGTGTATAAGAGACAG -b GTCTCGTGGGCTCGG -b AATGATACGGCGACCACCGAGATCTACAC -b CAAGCAGAAGACGGCATACGAGAT -b GACGCTGCCGACGA -b CTGTCTCTTATACACATCT -b CCGAGCCCACGAGAC -b GTGTAGATCTCGGTGGTCGCCGTATCATT -b ATCTCGTATGCCGTCTTCTGCTTG -o "../01c_cutadapt_merged_qualityfiltered/""$line""_BerryCrust.collapsed.AR" "$line""_BerryCrust.collapsed" ; done < ../Sample_list.txt

while read line; do    
cutadapt --times 10 -b TCGTCGGCAGCGTC -b AGATGTGTATAAGAGACAG -b GTCTCGTGGGCTCGG -b AATGATACGGCGACCACCGAGATCTACAC -b CAAGCAGAAGACGGCATACGAGAT -b GACGCTGCCGACGA -b CTGTCTCTTATACACATCT -b CCGAGCCCACGAGAC -b GTGTAGATCTCGGTGGTCGCCGTATCATT -b ATCTCGTATGCCGTCTTCTGCTTG -o "../01c_cutadapt_merged_qualityfiltered/""$line""_MiFish.collapsed.AR" "$line""_MiFish.collapsed" ; done < ../Sample_list.txt

cd ../01c_cutadapt_merged_qualityfiltered/
fastqc * -o ../01d_fastqc_merged_filtered_adapter_removed/
cd ../01d_fastqc_merged_filtered_adapter_removed
multiqc .

Yay!  This looks much better!  Now just NC45 Berry Crust and MiFish and NC 60 still have some traces of adapters still in them so I'll download them and look

OK, I can't actually find these adapters in these sequences like I could before, so perhaps some of the actual sequences are being considered adapter-like sequences by fastqc/multiqc
Done with adapter removal!


Now to try to put in the same annotations as ngsfilter does.  At a minimum we will need the sample name and the assay (MiFish or BerryCrust), ideally length too
from Daniel Marquina's eDNA tutorial: https://metagusano.github.io/publications/Bioinformatic%20Pipeline%20For%20Metabarcoding.pdf
(he did his immediately, before trimming or merging)
for file in *_R1.fastq; do
sample=${file#*dir/}
sample=${sample%_R1.fastq}
obiannotate -S sample:${sample} ${sample}_R1.fastq > ${sample}_R1.annotated.fastq
obiannotate -S sample:${sample} ${sample}_R2.fastq > ${sample}_R2.annotated.fastq 
done

Annotations:
for file in *_BerryCrust.collapsed.AR; do
sample=${file#*dir/}
sample=${sample%_BerryCrust.collapsed.AR}
obiannotate -S sample:${sample} -S assay:BerryCrust --length ${sample}_BerryCrust.collapsed.AR > ${sample}_BerryCrust.collapsed.AR.annotated.fastq
done

for file in *_MiFish.collapsed.AR; do
sample=${file#*dir/}
sample=${sample%_MiFish.collapsed.AR}
obiannotate -S sample:${sample} -S assay:MiFish --length ${sample}_MiFish.collapsed.AR > ${sample}_MiFish.collapsed.AR.annotated.fastq
done

ANNOTATIONS ADDED!

Now I'd like to filter for sequences between 140-260 bp long.  Also keep quality over 20, maximum 0 Ns and 1 expected errors.
mkdir ../01e_trimmed_QF
mkdir ../01e_trimmed_QF/discards
while read line; do
cutadapt -q 20 --max-n 0 --max-ee 1 -m 140 -M 260 --too-short-output "../01e_trimmed_QF/discards/""$line""_BerryCrust.under140bp" --too-long-output "../01e_trimmed_QF/discards/""$line""_BerryCrust.over260bp" -o "../01e_trimmed_QF/""$line""_BerryCrust.final.annotated.fastq" "$line""_BerryCrust.collapsed.AR.annotated.fastq"; done < ../Sample_list.txt
while read line; do
cutadapt -q 20 --max-n 0 --max-ee 1 -m 140 -M 260 --too-short-output "../01e_trimmed_QF/discards/""$line""_MiFish.under140bp" --too-long-output "../01e_trimmed_QF/discards/""$line""_MiFish.over260bp" -o "../01e_trimmed_QF/""$line""_MiFish.final.annotated.fastq" "$line""_MiFish.collapsed.AR.annotated.fastq"; done < ../Sample_list.txt

cd ../01e_trimmed_QF/
fastqc * -o ../01f_fastqc_trimmed_QF_final_annotated/
cd ../01f_fastqc_trimmed_QF_final_annotated
multiqc .

Looks great!  they're ready for the next steps now


OK, so according to the Physalia metabarcoding pipeline, it's time to combine them all and count how many sequences are in each for some reason.  
I'll do this for the MiFish and BerryCrust together, in 02_combine_and_count_all but I can always go back and do them separately if this doesn't work.  
But not yet.  I'll keep following the eDNAFlow pipeline and see what happens even though I don't see a chimera removal step.

According to the eDNAFlow pipeline, it's time to start Process 5: Relabel file for usearch.
Their input for example, is called F18_277.fastq and the output is: F18_277.relabeled.fastq  But the headers are gone and all of the associated metadata that was there before.  Now it's just 18_277.1 18_277.2 etc in the fastq sample names.
script:
   if(mode == 'usearch32')
	   """
	   for files in ${fastqs}
	   do
	   label=\$(echo \$files | cut -d '/' -f 3 | cut -d '.' -f 1)
	   usearch -fastx_relabel \$files -prefix \${label}. -fastqout \${label}.relabeled.fastq 
	   done
	 
	   for files in *.relabeled.fastq
	   do
	   name=\$(echo \$files | cut -d '/' -f '2' | cut -d '.' -f 1)
	   echo \${name} >> CountOfSeq.txt
	   grep "^@\${name}" \$files | wc -l >> CountOfSeq.txt
	   done 
	   cat *.relabeled.fastq > "${sample_id}_QF_Dmux_minLF_relabeled4Usearch.fastq"
	  
	   usearch -fastx_get_sample_names *_relabeled4Usearch.fastq -output sample.txt
	   usearch -fastq_filter *_relabeled4Usearch.fastq -fastaout ${sample_id}.fasta
	   
	   awk '/^>/ {print(\$0)}; /^[^>]/ {print(toupper(\$0))}' *.fasta > ${sample_id}_upper.fasta
	 
	   """
first, to test this, I will do the label part and view what it says.  

for files in *
do
label=$(echo ${files} | cut -d '/' -f '2' | cut -d '.' -f 1)
echo ${label}
done

This code worked.  Ok, the labels are both the sample and assay, just everything before the . so that's perfect.  It needed to be re-written a bit to work.

downloaded usearch32 11.0.667 https://www.drive5.com/usearch/download.html gunzipped and made it excecutable with chmod 777 and renamed it like this: mv usearch11.0.667_i86linux32 usearch
in eldridge@leopardus:/storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data

for files in *
do
label=$(echo ${files} | cut -d '/' -f '2' | cut -d '.' -f 1)
../usearch -fastx_relabel ${files} -prefix ${label}. -fastqout ../05_relabel_for_usearch/${label}.relabeled.fastq -keep_annots
done
 
this worked when run from 01e_trimmed_QF directory, but the sample= annotation was removed. even though the original code removed all of the annotations, this makes me nervous.  
I want to add sample back in.  Even though this may be unnecessary.

#run from 05_relabel_for_usearch
for file in *_MiFish.relabeled.fastq; do
sample=${file#*dir/}
sample=${sample%_MiFish.relabeled.fastq}
obiannotate -S sample:${sample} ${sample}_MiFish.relabeled.fastq > ${sample}_MiFish.relabeled_annotated.fastq
done

for file in *_BerryCrust.relabeled.fastq; do
sample=${file#*dir/}
sample=${sample%_BerryCrust.relabeled.fastq}
obiannotate -S sample:${sample} ${sample}_BerryCrust.relabeled.fastq > ${sample}_BerryCrust.relabeled_annotated.fastq
done

OK, next part of the relabeling!  
It looks like this concatenates everything back into the original "sample name" which is really the library name.
Then it counts how many reads are in each actual sample and puts it in CountofSeq.txt
I'm going to keep the MiFish and BerryCrust separate so they can each be processed with with their respective databases.

Not working with their code, so switching to my way:

while read line; do  
echo "$line""_MiFish" >> CountOfSeq.txt
grep "@""$line" "$line""_MiFish.relabeled_annotated.fastq" | wc -l >> CountOfSeq.txt; done < ../Sample_list.txt	   
while read line; do  
echo "$line""_BerryCrust" >> CountOfSeq.txt
grep "@""$line" "$line""_BerryCrust.relabeled_annotated.fastq" | wc -l >> CountOfSeq.txt; done < ../Sample_list.txt

while read line; do  
echo "$line""_BerryCrust" >> BerryCrust_CountOfSeq.txt
grep "@""$line" "$line""_BerryCrust.relabeled_annotated.fastq" | wc -l >> BerryCrust_CountOfSeq.txt; done < ../Sample_list.txt

while read line; do  
echo "$line""_MiFish" >> MiFish_CountOfSeq.txt
grep "@""$line" "$line""_MiFish.relabeled_annotated.fastq" | wc -l >> MiFish_CountOfSeq.txt; done < ../Sample_list.txt

### And the final part of their code edited for my purpose: but keeping the MiFish and BerryCrust separate so they can each be processed with their separate databases.

cat *_MiFish.relabeled_annotated.fastq > NC-MiFish_annotated_relabeled4Usearch.fastq
	  
	   ../usearch -fastx_get_sample_names NC-MiFish_annotated_relabeled4Usearch.fastq -output MiFishsample.txt

cat *_BerryCrust.relabeled_annotated.fastq > NC-BerryCrust_annotated_relabeled4Usearch.fastq	   
	   ../usearch -fastx_get_sample_names NC-BerryCrust_annotated_relabeled4Usearch.fastq -output BerryCrustsample.txt
	   
cat NC-MiFish_annotated_relabeled4Usearch.fastq NC-BerryCrust_annotated_relabeled4Usearch.fastq > NC-combined_annotated_relabeled4Usearch.fastq
../usearch -fastx_get_sample_names NC-combined_annotated_relabeled4Usearch.fastq -output combined_NC_samples.txt



changed their usearch conversion to fasta and awk to uppercase command for:
conda activate genomics
seqtk seq -A -U NC-MiFish_annotated_relabeled4Usearch.fastq > NC-MiFish_annotated.fasta  
worked!  -A puts it in FASTA format, -N removes ALL Ns -U makes it uppercase  I already removed all Ns so I won't do that again.

seqtk seq -A -U NC-BerryCrust_annotated_relabeled4Usearch.fastq > NC-BerryCrust_annotated.fasta
seqtk seq -A -U NC-combined_annotated_relabeled4Usearch.fastq > NC-combined_annotated.fasta

I did the combined version just in case I want that later, and if annotations get in the way I can come back and re-do from line 248 of this file, leaving out the --keep_annots and the re-annotating the sample= field

Now it's time to drop in to eDNAFlow with demultiplexed data

First, I have to make sure my blastdb from crabs is formatted correctly.

See Word Document NC_Sharks_data_steps (or similar which is mostly about the CRABS database steps)

OK, I think I've got it.
/storage/eldridge/Database_for_fish_and_crustaceans/MiFish_blastdb

From the eDNAFlow documentation:
nextflow run eDNAFlow.nf --barcode 'pe_bc*' --blast_db 'Path2TestBlastDataset/file.fasta' --custom_db 'path2/customDatabase/myDb' [OPTIONS]
--skipDemux: It's a boolean

--demuxedInput 'demuxedFile.fasta': provide name of the fasta file holding all the demultiplexed sequences. Format of the sample identifier must match USEARCH requirements.

    At least one of the below databases must be specifieThe --blast_db 'absolutePath2/LocalGenbankDatabase/nt': the absolute path to where nt databse is stored

--custom_db 'absolutePath2/customDatabase/myDb': the absolute path to where custom database is stored

From my history playing around with installing and running the test_data, I already have the singularity images stored so it doesn't have to download them every time:
nextflow run eDNAFlow.nf --singularityDir "~/eDNAFlow/install/singularity_images/"

ok, let's look try it!


conda activate eDNAflow
put eDNAFlow scripts into the folder with the data... /storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/demuxed_eDNAFlow
"Make sure eDNAFlow scripts (including eDNAFlow.nf, nextflow.config and lulu.R), conf and LCA_taxonomyAssignment_scripts folders are in the same directory where your unzipped sequencing and Multiplex identifier (MID) tag (here defined as “barcode”) files exist."

-from the eDNAFlow documentation: https://github.com/mahsa-mousavi/eDNAFlow#basic-command-usage

mkdir demuxed_eDNAFlow
cp ../05_relabel_for_usearch/NC-MiFish_annotated.fasta .
cp -r ~/eDNAFlow .
nextflow run eDNAFlow/eDNAFlow.nf --skipDemux --demuxedInput NC-MiFish_annotated.fasta --custom_db /storage/eldridge/Database_for_fish_and_crustaceans/MiFish_blastdb

started!  but wants the blastdb stuff all in this folder too
cp /storage/eldridge/Database_for_fish_and_crustaceans/MiFish_blastdb* .
nextflow run eDNAFlow/eDNAFlow.nf --skipDemux --demuxedInput NC-MiFish_annotated.fasta --custom_db MiFish_blastdb --singularityDir "/storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/demuxed_eDNAFlow/work/singularity" --resume

there was something wrong with the blastdb 

trying agan with the re-made one:
nextflow run eDNAFlow/eDNAFlow.nf --skipDemux --demuxedInput NC-MiFish_annotated.fasta --custom_db /storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/demuxed_eDNAFlow/MiFish_blastdb

blast is still failing with my custom database.

here's the command from the work directory:
blastn -task blastn -db "/storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/demuxed_eDNAFlow/MiFish_blastdb" -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" -perc_identity 95 -evalue 1e-3 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -qcov_hsp_perc 100 -max_target_seqs 10 -query NC-MiFish_annotated.fasta_zotus.fasta -out NC-MiFish_annotated.fasta_blast_Result.tab -num_threads 1

it told me to download wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
when I was in the work folder troubleshooting this, so I did and tar -xvf taxdb.tar.gz got it unpacked.  
The above command worked!  but only in the work folder d2!

OK, I'll just continue following the eDNAFlow pipeline on my own since there are only 2 steps left!
06 finished and 07 is partially finished.

so i copied the contents of d2/72b5ad3ea274674b7ea4500b09e4bf to demuxed_eDNAFlow/07_blast_NC-MiFish_annotated.fasta/ and out of the demuxed_eDNAFlow directory where I was trying to make eDNAFlow work, and into one directory up where I have just been following along with eDNAFlow.
/storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/07_blast_NC-MiFish_annotated.fasta

I'm not sure all of 07 has been completed so I'm going to pick up at 06_unique_zotus_usearch
mkdir 06_unique_zotus_usearch
cd 06_unique_zotus_usearch

../usearch -fastx_uniques ../demuxed_eDNAFlow/NC-MiFish_annotated.fasta -sizeout -fastaout NCNC-MiFish_annotated_Unq.fasta
results:
Fish_annotated.fasta -sizeout -fastaout NCNC-MiFish_annotated_Unq.fasta
usearch v11.0.667_i86linux32, 4.0Gb RAM (132Gb total), 64 cores
(C) Copyright 2013-18 Robert C. Edgar, all rights reserved.
https://drive5.com/usearch

License: personal use only

00:04 1.0Gb   100.0% Reading ../demuxed_eDNAFlow/NC-MiFish_annotated.fasta
00:04 1.0Gb  CPU has 64 cores, defaulting to 10 threads                   
00:05 1.8Gb   100.0% DF
00:05 1.8Gb  3438227 seqs, 97931 uniques, 72297 singletons (73.8%)
00:05 1.8Gb  Min size 1, median 1, max 456872, avg 35.11
00:05 1.2Gb   100.0% Writing NCNC-MiFish_annotated_Unq.fasta

oops, too many NCs in the name:
mv NCNC-MiFish_annotated_Unq.fasta NC-MiFish_annotated_Unq.fasta

../usearch -unoise3 NC-MiFish_annotated_Unq.fasta  -zotus NC-MiFish_annotated_zotus.fasta -tabbedout NC-MiFish_annotated_Unq_unoise3.txt -minsize 8    
results:
00:01 69Mb    100.0% Reading NC-MiFish_annotated_Unq.fasta
00:01 38Mb      0.0% 0 amplicons, 0 bad (size >= 456872)  
WARNING: Shifted sequences detected

00:01 44Mb    100.0% 236 amplicons, 603375 bad (size >= 8)
00:01 51Mb    100.0% 81 good, 155 chimeras                
00:01 51Mb    100.0% Writing zotus

hmmm, not loving this result.  I want to do the chimera removal like the obitools pipeline suggests instead.  But I could also change the minsize
For now I'll keep it and see what it looks like in the end.  Not sure what shifted sequences mean.

continuing on...
../usearch -otutab ../demuxed_eDNAFlow/NC-MiFish_annotated.fasta -zotus NC-MiFish_annotated_zotus.fasta -otutabout zotuTable.txt -mapout zmap.txt
results:
00:00 41Mb    100.0% Reading NC-MiFish_annotated_zotus.fasta
00:00 7.1Mb   100.0% Masking (fastnucleo)                   
00:00 7.9Mb   100.0% Word stats          
00:00 7.9Mb   100.0% Alloc rows
00:00 7.9Mb   100.0% Build index
00:00 41Mb   CPU has 64 cores, defaulting to 10 threads
02:18 137Mb   100.0% Searching NC-MiFish_annotated.fasta, 99.7% matched 
3428925 / 3438227 mapped to OTUs (99.7%)                               
02:18 137Mb  Writing zotuTable.txt
02:18 137Mb  Writing zotuTable.txt ...done.

cd ..
mkdir 07_blast
cd 07_blast
cp ../demuxed_eDNAFlow/taxdb.* .

blastn -db "/storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/demuxed_eDNAFlow/MiFish_blastdb" -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" -perc_identity 95 -evalue 1e-3 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -qcov_hsp_perc 100 -max_target_seqs 10 -query ../06_unique_zotus_usearch/NC-MiFish_annotated_zotus.fasta -out NC-MiFish_annotated_blast_Result.tab -num_threads 1

worked!

makeblastdb -in ../06_unique_zotus_usearch/NC-MiFish_annotated_zotus.fasta -parse_seqids -dbtype nucl -out NC-MiFish_annotated_zotus
results:
Building a new DB, current time: 11/05/2022 10:49:05
New DB name:   /storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/07_blast/NC-MiFish_annotated_zotus
New DB title:  ../06_unique_zotus_usearch/NC-MiFish_annotated_zotus.fasta
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 3000000000B
Adding sequences from FASTA; added 81 sequences in 0.057132 seconds.

blastn -db NC-MiFish_annotated_zotus \
             -outfmt "6 qseqid sseqid pident" \
             -out match_list.txt -qcov_hsp_perc 80 \
             -perc_identity 84 -query ../06_unique_zotus_usearch/NC-MiFish_annotated_zotus.fasta \
             -num_threads 1
             
worked!  now I have a list of similar zotus 84% similarity threshold.

cd ..
mkdir 08_lulu
cd 08_lulu/

cd following_eDNAFlow_w_demuxed_data/demuxed_eDNAFlow/eDNAFlow
cp lulu.R ../../08_lulu/

copied the script out of lulu.R in the eDNAFlow folder and made a new script incase I have to change it to use it outside of Nextflow.  scp ../../following_eDNAFlow_method/eDNAFlowlulu.R eldridge@leopardus.snrenet.arizona.edu:/storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/08_lulu
chmod 777 eDNAFlowlulu.R
it already has the names of the output files from step 6 and 7 hardcoded, so I'll just copy them to 08
cp ../07_blast/match_list.txt .
cp ../06_unique_zotus_usearch/zotuTable.txt .

from the eDNAFlow commands:
Rscript $lulu ${minMatch_lulu} (84 is the default value they use)
Rscript eDNAFlowlulu.R 84
first I need to install lulu!

library(devtools)
install_github("tobiasgf/lulu")  

but first apparently I have to install devtools.

with these three lines added to the beginning of the script, it started downloading devtools
install.packages("devtools")
library(devtools)
install_github("tobiasgf/lulu")

I made a conda environment just in case it changes anything where I am.  I'm in conda env lulu.R
I could have just done this part on my own computer, but I'm interested in learning how to do it on leopardus
but it says: Error in library(devtools) : there is no package called ‘devtools’
Execution halted even though it just did a lot of downloading... maybe it needs quotes?

OK, I got tired of trying to make it work on leopardus and I did it in RStudio on my computer. it took 2 minutes
All input and output folders (RStudio working directory) are in:
/Users/Eldridge/Desktop/Savannah/following_eDNAFlow_method

then I put them back onto leopardus to continue with step 9 (LCA assignment script)
mkdir ../09_taxonomy_assignment_LCA
cd ../09_taxonomy_assignment_LCA
cp ../08_lulu/curated_zotuTable.tab .
cp ../07_blast/NC-MiFish_annotated_blast_Result.tab .
cp ../06_unique_zotus_usearch/zotuTable.txt .


from the help command:
nextflow run ../eDNAFlow.nf --help
For running LCA taxonomy assignment script:
nextflow run eDNAFlow.nf --taxonomyAssignment --zotuTable "path2/curatedOruncurated_ZotuTable_file" --blastFile "path2/blastResult_file" --lca_output "my_lca_result" [OPTIONS]

nextflow run /storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/demuxed_eDNAFlow/eDNAFlow/eDNAFlow.nf --taxonomyAssignment --zotuTable zotuTable.txt --blastFile NC-MiFish_annotated_blast_Result.tab --lca_output "NC-MiFish_crabsdb_nolulu_results"
the first time I tried this, it failed because it couldn't find the scripts in this folder so I copied them in: 
cp -r ../demuxed_eDNAFlow/eDNAFlow/LCA_taxonomyAssignment_scripts/ .

ended with this error:
IndexError: list index out of range

This will take some more work to figure it out.  I'll have to go into the python scripts.

What I think is going on is that the tables I have don't have all of the fields they're supposed to.  

I should re-do step 6 and 7 to make sure they're right.  
------
OR, I could try it with the full nt blast database, but leopardus can't hold that much.  Loaded it on my new mac and brought over the results of step 6 to do step 7 blast on
my computer doesn't like usearch so I can't re-do step 6 on here.  I could on leopardus... maybe something to do with the new M1 max chip in the new mac.
and finally, it could be that I kept the annotations, and I should try to re-do it with the CRABs database on leopardus without the annotations.

So.  I first want to try blasting both MiFish and BerryCrust together against the entire nt database on my new mac.  But I need to do the usearch step on leopardus first.  
So for now I'll just do the MiFish.
Here's the command that worked (I think) on leopardus:
blastn -db "/storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/demuxed_eDNAFlow/MiFish_blastdb" -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" -perc_identity 95 -evalue 1e-3 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -qcov_hsp_perc 100 -max_target_seqs 10 -query ../06_unique_zotus_usearch/NC-MiFish_annotated_zotus.fasta -out NC-MiFish_annotated_blast_Result.tab -num_threads 1
now edited for my new mac:
blastn -db "/Users/Eldridge/Desktop/Savannah/following_eDNAFlow_method/blast_db/nt" -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" -perc_identity 95 -evalue 1e-3 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -qcov_hsp_perc 100 -max_target_seqs 10 -query ../06_unique_zotus_usearch/NC-MiFish_annotated_zotus.fasta -out NC-MiFish_annotated_blast_Result_2.tab -num_threads 1
results:
BLAST Database error: Error pre-fetching sequence data 

OK, not working.  
Back to leopardus:  I'll do the not-annotated MiFish against the CRABS database, and see if the annotations were the problem.

mkdir 05b_relabel_for_usearch_no_annotations
cd ../01e_trimmed_QF

for files in *
do
label=$(echo ${files} | cut -d '/' -f '2' | cut -d '.' -f 1)
../usearch -fastx_relabel ${files} -prefix ${label}. -fastqout ../05b_relabel_for_usearch_no_annotations/${label}.relabeled.fastq
done

mkdir summaryfiles

while read line; do  
echo "$line""_BerryCrust" >> summaryfiles/BerryCrust_CountOfSeq.txt
grep "@""$line" "$line""_BerryCrust.relabeled.fastq" | wc -l >> summaryfiles/BerryCrust_CountOfSeq.txt; done < ../Sample_list.txt

while read line; do  
echo "$line""_MiFish" >> summaryfiles/MiFish_CountOfSeq.txt
grep "@""$line" "$line""_MiFish.relabeled.fastq" | wc -l >> summaryfiles/MiFish_CountOfSeq.txt; done < ../Sample_list.txt

### And the final part of their code edited for my purpose: but keeping the MiFish and BerryCrust separate so they can each be processed with their separate databases.

cat *_MiFish.relabeled.fastq > summaryfiles/NC-MiFish_relabeled4Usearch.fastq
	  
	   ../usearch -fastx_get_sample_names summaryfiles/NC-MiFish_relabeled4Usearch.fastq -output summaryfiles/MiFishsamples.txt

cat *_BerryCrust.relabeled.fastq > summaryfiles/NC-BerryCrust_relabeled4Usearch.fastq	   
	   ../usearch -fastx_get_sample_names summaryfiles/NC-BerryCrust_relabeled4Usearch.fastq -output summaryfiles/BerryCrustsamples.txt

cd summaryfiles/	   
cat NC-MiFish_relabeled4Usearch.fastq NC-BerryCrust_relabeled4Usearch.fastq > NC-combined_relabeled4Usearch.fastq
../../usearch -fastx_get_sample_names NC-combined_relabeled4Usearch.fastq -output combined_NC_samples.txt



changed their usearch conversion to fasta and awk to uppercase command for:
conda activate genomics
seqtk seq -A -U NC-MiFish_relabeled4Usearch.fastq > NC-MiFish.fasta  
worked!  -A puts it in FASTA format, -N removes ALL Ns -U makes it uppercase  I already removed all Ns so I won't do that again.

seqtk seq -A -U NC-BerryCrust_relabeled4Usearch.fastq > NC-BerryCrust.fasta

seqtk seq -A -U NC-combined_relabeled4Usearch.fastq > NC-combined.fasta

eDNAFlow attempt with non_annotated files:
mkdir demuxed_eDNAFlow_no_annotations
cd demuxed_eDNAFlow
cp -r eDNAFlow/ ../demuxed_eDNAFlow_no_annotations/
cp FINAL_MiFish_NOTAX.fasta ../demuxed_eDNAFlow_no_annotations/
cp taxdb* ../demuxed_eDNAFlow_no_annotations/
cd ../demuxed_eDNAFlow_no_annotations

makeblastdb -in FINAL_MiFish_NOTAX.fasta -parse_seqids -dbtype nucl -taxid_map ../../../../Database_for_fish_and_crustaceans/nucl_gb.accession2taxid -out MiFish_blastdb -title MiFish_blastdb
results:
Building a new DB, current time: 11/09/2022 13:43:00
New DB name:   /storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/demuxed_eDNAFlow_no_annotations/MiFish_blastdb
New DB title:  MiFish_blastdb
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 3000000000B
Adding sequences from FASTA; added 10428 sequences in 0.192292 seconds.

I'll try eDNA_Flow but then I can always do it by hand.  Also I made note earlier that I want to try chimaera removal as well, so I need to try that method too.
nextflow run eDNAFlow/eDNAFlow.nf --skipDemux --demuxedInput NC-MiFish.fasta --custom_db MiFish_blastdb --singularityDir ~/eDNAFlow/install/singularity_images/

blast step still not working with eDNAFlow even though I removed the space that I thought might have been the issue in the Nextflow.nf file
mkdir 07_blast (in demuxed_eDNAFlow_no_annotations directory)
blastn -db "/storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/demuxed_eDNAFlow_no_annotations/MiFish_blastdb" -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" -perc_identity 95 -evalue 1e-3 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -qcov_hsp_perc 100 -max_target_seqs 10 -query 06_Uniques_ZOTUs_NC-MiFish.fasta/NC-MiFish.fasta_zotus.fasta -out 07_blast/NC-MiFish_blast_Result.tab -num_threads 1
worked very quickly but all of the taxonomy fields are still N/A

trying to make the database over again:
makeblastdb -in ../FINAL_MiFish_NOTAX.fasta -parse_seqids -dbtype nucl -taxid_map ../../../../../Database_for_fish_and_crustaceans/nucl_gb.accession2taxid -out /storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/demuxed_eDNAFlow_no_annotations/blastdb/MiFish_blastdb -title "MiFish_blastdb"


Building a new DB, current time: 11/09/2022 16:56:31
New DB name:   /storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/demuxed_eDNAFlow_no_annotations/blastdb/MiFish_blastdb
New DB title:  MiFish_blastdb
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 3000000000B
Adding sequences from FASTA; added 10428 sequences in 0.252965 seconds.

export BLASTDB=/storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/demuxed_eDNAFlow_no_annotations/blastdb
ls $BLASTDB/taxdb.*
blastn -task blastn -db /storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/demuxed_eDNAFlow_no_annotations/blastdb/MiFish_blastdb -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" -perc_identity 95 -evalue 1e-3 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -qcov_hsp_perc 100 -max_target_seqs 10 -query NC-MiFish.fasta_zotus.fasta -out 07_blast/NC-MiFish_blast_Result.tab -num_threads 1

results: all taxonomy is NA.  there are blastIDs though.  So I could try using that program taht will convert them and then figure out a way to insert them into the table properly...

OK, screw this, I'm going to do it a different way (more similar to the Physalia courses Owen Wangensteen way)

ORRRR dada2 in R because crabs has two output options for dada2
ok dada2 installed on my new macbook
following this tutorial: https://benjjneb.github.io/dada2/tutorial.html

worked!  but only 17 sequences assigned to species (ASV and dada2 taxonomy assignment, requires 100% match for species -more suitable algorithm for bacteria)

ok, OBItools2 method:
starting from 01e_trimmed_QF

cd /Users/Eldridge/Desktop/Savannah/NC_Sharks_obitools/2nd_try_demuxed_cleaned
scp -r eldridge@leopardus.snrenet.arizona.edu:/storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/01e_trimmed_QF/ .
mkdir 02_obiconvert_to_fasta

while read line; do  
obiconvert --sanger --fasta-output \
"$line""_MiFish.final.annotated.fastq" > "../02_obiconvert_to_fasta/""$line""_MiFish""".fasta""; done < ../Sample_list.txt	

while read line; do  
obiconvert --sanger --fasta-output \
"$line""_BerryCrust.final.annotated.fastq" > "../02_obiconvert_to_fasta/""$line""_BerryCrust""".fasta""; done < ../Sample_list.txt


cat *_MiFish.fasta  > MiFish.allsamples.fasta
obistat -c sample -a seq_length MiFish.allsamples.fasta > sample_stats_MiFish.txt 
obisplit -t sample MiFish.allsamples.fasta 
while read line; do  
rm "$line"".fasta"; done < ../Sample_list.txt


cat *_BerryCrust.fasta  > BerryCrust.allsamples.fasta
obistat -c sample -a seq_length BerryCrust.allsamples.fasta > sample_stats_BerryCrust.txt 

cd ..
mkdir owitools
cd owitools
cp -r ~/Desktop/Metabarcoding_Physalia/ubuntu/owi_tools/ .
chmod 777 *

for i in *_MiFish.fasta ; do obiuniq $i >  ${i/.fasta/.unique.fasta} & done
which RScript 
result on my computer: /usr/local/bin/RScript
copy and paste that into the first line of the owi_script that I want to use
also had to go into RStudio and install.packages("optparse")
for i in *_MiFish.unique.fasta ; do owi_obisample2vsearch -i $i & done
for i in *_MiFish.unique.vsearch.fasta ; do vsearch --uchime_denovo $i --sizeout --minh 0.90 --nonchimeras ${i/.unique.vsearch.fasta/.nonchimeras.fasta} --chimeras ${i/.unique.vsearch.fasta/.chimeras.fasta} --uchimeout ${i/.unique.vsearch.fasta/.uchimeout.txt} & done

#Concatenating non chimaeras in a single file
cat  *_MiFish.nonchimeras.fasta > MiFish.nonchimeras.fasta

#Returning to obitools format
../owitools/owi_vsearch2obifasta -i MiFish.nonchimeras.fasta
#Dereplicating sequences in MiFish.nonchimeras.fasta file
obiuniq -m sample MiFish.nonchimeras.vsearch.fasta > MiFish.unique.fasta

#Change identifiers by a short index
obiannotate --seq-rank MiFish.unique.fasta | obiannotate --set-identifier '"MiFish_%09d" % seq_rank' > MiFish.new.fasta

#obtain a table of abundances
obitab -o MiFish.new.fasta >  MiFish.new.tab
#end of hands_on_session3
#back to vsearch format
../owitools/owi_obifasta2vsearch -i MiFish.new.fasta -o MiFish.vsearch.fasta
#cluster with swarm
conda activate eDNA
swarm -d 4 -z -t 40 -o MiFish.SWARM4nc_output -s MiFish.SWARM4nc_stats -w MiFish.SWARM4nc_seeds.fasta MiFish.vsearch.fasta

error duplicated sequence identifier: MiFish_0000000

#I could change the number of leading zeros since it only seems to recognize the first 7
obiannotate --seq-rank MiFish.unique.fasta | obiannotate --set-identifier '"MiFish_%07d" % seq_rank' > MiFish.new.fasta

#obtain a table of abundances
obitab -o MiFish.new.fasta >  MiFish.new.tab
#end of hands_on_session3
#back to vsearch format
../owitools/owi_obifasta2vsearch -i MiFish.new.fasta -o MiFish.vsearch.fasta
#cluster with swarm
conda activate eDNA
swarm -d 4 -z -t 40 -o MiFish.SWARM4nc_output -s MiFish.SWARM4nc_stats -w MiFish.SWARM4nc_seeds.fasta MiFish.vsearch.fasta
for i in 2 3 5 6 7 8 9 10 11 12 13 14; do
swarm -d $i -z -t 40 -o MiFish.SWARM$inc_output -s MiFish.SWARM$inc_stats -w MiFish.SWARM$inc_seeds.fasta MiFish.vsearch.fasta; done
swarm -d 1 -f -z -t 40 -o MiFish.SWARM1nc_output -s MiFish.SWARM1nc_stats -w MiFish.SWARM1nc_seeds.fasta MiFish.vsearch.fasta
it appears that the only stats files are the ones I didn't run in the loop
swarm -d 5 -z -t 40 -o MiFish.SWARM5nc_output -s MiFish.SWARM5nc_stats -w MiFish.SWARM5nc_seeds.fasta MiFish.vsearch.fasta; done

#I plotted the number of swarms vs d values and it starts to even out at d=5 (713 swarms) so that's what I'll use first

../owitools/owi_recount_swarm MiFish.SWARM5nc_output MiFish.new.tab
#Remove singletons
#sed -i 's/;size/ size/g' MiFish.SWARM5nc_seeds.fasta #this doesn't work so I'm going to use owi tools to go from vsearch format back to obitools format
 ../owitools/owi_vsearch2obifasta -i MiFish.SWARM5nc_seeds.fasta -o MiFish.SWARM5nc_seeds.obi.fasta
obigrep -p 'count>1' MiFish.SWARM5nc_seeds.fasta > MiFish.seeds_nonsingleton.fasta

#now putting it on leopardus to analyze against the obi3 MiFish database I made there (see obitools3_commands.txt)
scp MiFish.seeds_nonsingleton.fasta eldridge@leopardus.snrenet.arizona.edu:/storage/eldridge/Database_for_fish_and_crustaceans/obitools

#on leopardus in obitools3
conda activate obi3
source obi3-env/bin/activate
#make reference database of vertebrates that amplify with these primers following obitools3 tutorial: https://git.metabarcoding.org/obitools/obitools3/-/wikis/Wolf-tutorial-with-the-OBITools3
mkdir EMBL
 cd EMBL

wget -nH --cut-dirs=6 -A 'STD_VRT*.dat.gz' -R 'STD_HUM*.dat.gz','STD_ENV*.dat.gz' -m -np ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/std/
#this only downloaded 4 files, although it looks like there are 17 in the ftp site
cd ..
obi import --embl EMBL fish/Miya_refs
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

obi import --taxdump taxdump.tar.gz fish/taxonomy/my_tax


obi ecopcr -e 3 -l 100 -L 260 -F GTCGGTAAAACTCGTGCCAGC -R CATAGTGGGGTATCTAATCCCAGTTTG --taxonomy fish/taxonomy/my_tax fish/Miya_refs fish/Miya_ecopcr_refs
obi grep --require-rank=species --require-rank=genus --require-rank=family --taxonomy fish/taxonomy/my_tax fish/Miya_ecopcr_refs fish/Miya_refs_clean

obi build_ref_db -t 0.97 --taxonomy fish/taxonomy/my_tax fish/Miya_refs_clean fish/Miya_refs_clean_97

 

obi import --nuc MiFish.seeds_nonsingleton.fasta fish/MiFish_swarm5_nucimport

obi ecotag -m 0.97 --taxonomy fish/taxonomy/my_tax -R fish/Miya_refs_clean_97 fish/MiFish_swarm5_nucimport fish/MiFish_swarm5_nucimp_tax_assigned

obi stats -c SCIENTIFIC_NAME fish/MiFish_swarm5_nucimp_tax_assigned

IT WORKED!
swarm d=5

SCIENTIFIC_NAME                 count   total
None                            207     1569968 
Sphyrna tiburo                  1       534336  
Carcharhinus acronotus          1       457646  
Leiostomus xanthurus            1       339782  
Brevoortia                      1       235799  
Carcharhinus                    1       53850   
Cynoscion regalis               1       43142   
Lagodon rhomboides              1       40346   
Peprilus                        1       27054   
Cynoscion nothus                1       20640   
Symphurus                       1       16119   
Synodus foetens                 1       15738   
Ophichthus gomesii              1       15517   
Myrophis punctatus              1       14285   
Xyrichtys novacula              1       13217   
Chilomycterus schoepfii         1       6734
Sciaenops ocellatus             1       5994
Trichiurus lepturus             1       5742
Hippocampus erectus             1       4567
Stephanolepis                   1       3524
Syngnathus louisianae           1       3229
Paralichthys lethostigma        1       1936
Pleuronectidae                  1       991
Trachinotus                     1       257
Pomatomus saltatrix             1       65
Oncorhynchus                    1       20
Chaetodipterus faber            1       10
Echeneis                        1       8
Lutjanus guttatus               1       7
Thunnus                         1       7
Scomberomorus maculatus         1       7
Sphoeroides parvus              1       6

#now trying other d-values 
mv MiFish.seeds_nonsingleton.fasta MiFish.5seeds_nonsingleton.fasta

swarm d=1
../owitools/owi_recount_swarm MiFish.SWARM1nc_output MiFish.new.tab
../owitools/owi_vsearch2obifasta -i MiFish.SWARM1nc_seeds.fasta -o MiFish.SWARM1nc_seeds.obi.fasta
obigrep -p 'count>1' MiFish.SWARM1nc_seeds.obi.fasta > MiFish.1seeds_nonsingleton.fasta

scp MiFish.1seeds_nonsingleton.fasta eldridge@leopardus.snrenet.arizona.edu:/storage/eldridge/Database_for_fish_and_crustaceans/obitools
obi import --nuc MiFish.1seeds_nonsingleton.fasta fish/MiFish_swarm1_nucimport
obi ecotag -m 0.97 --taxonomy fish/taxonomy/my_tax -R fish/Miya_refs_clean_97 fish/MiFish_swarm1_nucimport fish/MiFish_swarm1_nucimp_tax_assigned
obi stats -c SCIENTIFIC_NAME fish/MiFish_swarm1_nucimp_tax_assigned
SCIENTIFIC_NAME                 count   total
None                            1138    1560354 
Sphyrna tiburo                  7       530918  
Carcharhinus acronotus          10      454628  
Leiostomus xanthurus            8       339814  
Brevoortia                      10      234214  
Carcharhinus                    10      61269   
Cynoscion regalis               5       45840   
Lagodon rhomboides              4       39890   
Peprilus                        5       27023   
Cynoscion nothus                6       24527   
Symphurus                       4       15944   
Synodus foetens                 3       15711   
Ophichthus gomesii              4       15468   
Myrophis punctatus              4       14249   
Xyrichtys novacula              4       13143   
Chilomycterus schoepfii         3       6728
Sciaenops ocellatus             4       6593
Trichiurus lepturus             4       5653
Hippocampus erectus             4       4543
Stephanolepis                   5       3511
Syngnathus louisianae           4       3190
Paralichthys lethostigma        3       1936
Pleuronectidae                  2       986
Trachinotus                     1       255
Pomatomus saltatrix             1       65
Sardinella aurita               3       35
Oncorhynchus                    1       20
Mugil thoburni                  2       15*
Chaetodipterus faber            1       10
Echeneis                        1       8
Lutjanus guttatus               1       7
Thunnus                         1       7
Sphoeroides parvus              1       6
Cynoscion nebulosus             1       6*
Scomberomorus maculatus         1       6*

#got 3 more species-level assignments I marked them with *
Cynoscion nebulosus and Scomberomorus maculatus occur in the correct area, but Mugil thoburni is from the pacific and galapagos

#now trying with swarm5 but no singleton removal
scp MiFish.SWARM5nc_seeds.fasta eldridge@leopardus.snrenet.arizona.edu:/storage/eldridge/Database_for_fish_and_crustaceans/obitools
obi import --nuc MiFish.SWARM5nc_seeds.fasta fish/MiFish_swarm5_wsingles
obi ecotag -m 0.97 --taxonomy fish/taxonomy/my_tax -R fish/Miya_refs_clean_97 fish/MiFish_swarm5_wsingles fish/MiFish_swarm5_wsingles_tax_assigned
obi stats -c SCIENTIFIC_NAME fish/MiFish_swarm5_wsingles_tax_assigned

SCIENTIFIC_NAME                 count   total
None                            680     680
Sphyrna tiburo                  1       1
Carcharhinus acronotus          1       1
Leiostomus xanthurus            1       1
Brevoortia                      1       1
Carcharhinus                    1       1
Cynoscion regalis               1       1
Lagodon rhomboides              1       1
Peprilus                        1       1
Cynoscion nothus                1       1
Symphurus                       1       1
Synodus foetens                 1       1
Ophichthus gomesii              1       1
Myrophis punctatus              1       1
Xyrichtys novacula              1       1
Chilomycterus schoepfii         1       1
Sciaenops ocellatus             1       1
Trichiurus lepturus             1       1
Hippocampus erectus             1       1
Stephanolepis                   1       1
Syngnathus louisianae           1       1
Paralichthys lethostigma        1       1
Pleuronectidae                  1       1
Trachinotus                     1       1
Pomatomus saltatrix             1       1
Oncorhynchus                    1       1
Chaetodipterus faber            1       1
Echeneis                        1       1
Lutjanus guttatus               1       1
Thunnus                         1       1
Scomberomorus maculatus         1       1
Sphoeroides parvus              1       1
Centropristis striata           1       1
Eucinostomus gula               1       1

# same order, two more viable species at the bottom of the list, weird numbering though: let's see if this fixes the numbers
obigrep -p 'count>0' MiFish.SWARM5nc_seeds.obi.fasta > MiFish.5seeds_wsingles.fasta
scp MiFish.5seeds_wsingles.fasta eldridge@leopardus.snrenet.arizona.edu:/storage/eldridge/Database_for_fish_and_crustaceans/obitools
obi import --nuc MiFish.5seeds_wsingles.fasta fish/MiFish_swarm5_wsingles_count
obi ecotag -m 0.97 --taxonomy fish/taxonomy/my_tax -R fish/Miya_refs_clean_97 fish/MiFish_swarm5_wsingles_count fish/MiFish_swarm5_wsingles_count_tax_assigned
obi stats -c SCIENTIFIC_NAME fish/MiFish_swarm5_wsingles_count_tax_assigned

SCIENTIFIC_NAME                 count   total
None                            680     1570441 
Sphyrna tiburo                  1       534336  
Carcharhinus acronotus          1       457646  
Leiostomus xanthurus            1       339782  
Brevoortia                      1       235799  
Carcharhinus                    1       53850   
Cynoscion regalis               1       43142   
Lagodon rhomboides              1       40346   
Peprilus                        1       27054   
Cynoscion nothus                1       20640   
Symphurus                       1       16119   
Synodus foetens                 1       15738   
Ophichthus gomesii              1       15517   
Myrophis punctatus              1       14285   
Xyrichtys novacula              1       13217   
Chilomycterus schoepfii         1       6734
Sciaenops ocellatus             1       5994
Trichiurus lepturus             1       5742
Hippocampus erectus             1       4567
Stephanolepis                   1       3524
Syngnathus louisianae           1       3229
Paralichthys lethostigma        1       1936
Pleuronectidae                  1       991
Trachinotus                     1       257
Pomatomus saltatrix             1       65
Oncorhynchus                    1       20
Chaetodipterus faber            1       10
Echeneis                        1       8
Lutjanus guttatus               1       7
Thunnus                         1       7
Scomberomorus maculatus         1       7
Sphoeroides parvus              1       6
Centropristis striata           1       1*
Eucinostomus gula               1       1*

#this added Centropristis striata and Eucinostomus gula which both occur in the correct region
may as well look at d=4 while we're at it:
swarm d=4
../owitools/owi_recount_swarm MiFish.SWARM1nc_output MiFish.new.tab
../owitools/owi_vsearch2obifasta -i MiFish.SWARM4nc_seeds.fasta -o MiFish.SWARM4nc_seeds.obi.fasta
obigrep -p 'count>1' MiFish.SWARM4nc_seeds.obi.fasta > MiFish.4seeds_nonsingleton.fasta

scp MiFish.4seeds_nonsingleton.fasta eldridge@leopardus.snrenet.arizona.edu:/storage/eldridge/Database_for_fish_and_crustaceans/obitools
obi import --nuc MiFish.4seeds_nonsingleton.fasta fish/MiFish_swarm4_nosingles
obi ecotag -m 0.97 --taxonomy fish/taxonomy/my_tax -R fish/Miya_refs_clean_97 fish/MiFish_swarm4_nosingles fish/MiFish_swarm4_nosingles_tax_assigned
obi stats -c SCIENTIFIC_NAME fish/MiFish_swarm4_nosingles_tax_assigned

SCIENTIFIC_NAME                 count   total
None                            250     1558082
Sphyrna tiburo                  1       534295
Carcharhinus acronotus          1       460028
Leiostomus xanthurus            1       339807
Brevoortia                      1       235797
Carcharhinus                    2       60265
Cynoscion regalis               1       44243
Lagodon rhomboides              1       40345
Peprilus                        1       27054
Cynoscion nothus                1       22094
Symphurus                       1       16009
Synodus foetens                 1       15738
Ophichthus gomesii              1       15515
Myrophis punctatus              1       14283
Xyrichtys novacula              1       13216
Chilomycterus schoepfii         1       6733
Sciaenops ocellatus             1       6444
Trichiurus lepturus             1       5741
Hippocampus erectus             1       4567
Stephanolepis                   1       3524
Syngnathus louisianae           1       3229
Paralichthys lethostigma        1       1936
Pleuronectidae                  1       991
Trachinotus                     1       257
Pomatomus saltatrix             1       65
Oncorhynchus                    1       20
Chaetodipterus faber            1       10
Echeneis                        1       8
Lutjanus guttatus               1       7
Thunnus                         1       7
Scomberomorus maculatus         1       7
Sphoeroides parvus              1       6

obi align -t 0.95 fish/MiFish_swarm4_nosingles_tax_assigned fish/aligned_swarm4_assigned_sequences  
            
obi export --fasta-output fish/MiFish_swarm4_nosingles_tax_assigned -o NCSharks_MiFish_results.fasta

obi export --tab-output fish/aligned_swarm4_assigned_sequences > NCSharks_MiFish_results.csv

cd /Users/Eldridge/Desktop/Savannah/NC_Sharks_obitools/results
scp eldridge@leopardus.snrenet.arizona.edu:/storage/eldridge/Database_for_fish_and_crustaceans/obitools/NCSharks_MiFish_results* .

code from this tutorial (but they did it just after swarm clustering): https://metagusano.github.io/publications/Bioinformatic%20Pipeline%20For%20Metabarcoding.pdf
vsearch --usearch_global NCSharks_MiFish_results.fasta --db NCSharks_MiFish_results.fasta --self --id .84 --iddef 1 --userout NCSharks_MiFish_results.match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
seemed to work!

I'll do the swarm4.seeds too like in this tutorial in case that works differently.
vsearch --usearch_global ../2nd_try_demuxed_cleaned/02_obiconvert_to_fasta/MiFish.SWARM4nc_seeds.obi.fasta --db ../2nd_try_demuxed_cleaned/02_obiconvert_to_fasta/MiFish.SWARM4nc_seeds.obi.fasta --self --id .84 --iddef 1 --userout MiFish.SWARM4nc_seeds.obi.match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
#owi_combine must be done before LULU
first,_recount_swarm!
../owitools/owi_recount_swarm MiFish.SWARM4nc_output MiFish.new.tab
#results:MiFish.SWARM4nc_output.counts.csv

#owi_combine -i ULOY.ecotag.fasta.annotated.csv -a ULOY.SWARM13nc_output.counts.csv -o ULOY.All_MOTUs.csv
owi_combine -i NCSharks_MiFish_results.csv -a MiFish.SWARM4nc_output.counts.csv -o FISH.All_MOTUs.csv

it's not liking the output format of obi3.

onto the swarm output LULU then obi3 taxonomic assignment
I decided to use swarm 4 because it resulted in the same number of motus being assigned to only 1 species as 5 did.
cp ../02_obiconvert_to_fasta/MiFish.4seeds_nonsingleton.fasta .
sed -i -e 's/ size=\([0-9]\+\);//g' MiFish.4seeds_nonsingleton.fasta
vsearch --usearch_global ../2nd_try_demuxed_cleaned/02_obiconvert_to_fasta/MiFish.SWARM4nc_seeds.obi.fasta --db ../2nd_try_demuxed_cleaned/02_obiconvert_to_fasta/MiFish.SWARM4nc_seeds.obi.fasta --self --id .84 --iddef 1 --userout MiFish.SWARM4nc_seeds.obi.match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10


#with this command we are transforming the output.counts.csv from the owi_recount step into a simpler table with only id and the reads/sample, which is the only info that LULU considers. The argument -f 1, A-B, indicates to print only the columns 1 (id) and A-B (samples). You have to substitute A and B for the first and last sample column of your table: if you have 10 samples starting in column 15, then A=15, B=24.

####### A and B must change to the number of the column containing your first sample and last sample respectively#####
cut -d ";" -f1,A-B output_of_recount_swarm.csv > name.counts_LULU.txt && sed -e 's/;/\t/g' name.counts_LULU.txt

cp ../02_obiconvert_to_fasta/MiFish.SWARM4nc_output.counts.csv .
#opened in excel and inserted a row to count the columns.a=4,B=160
cut -d ";" -f1,4-160 MiFish.SWARM4nc_output.counts.csv > MiFish.SWARM4.counts_LULU.txt && sed -i -e 's/;/\t/g' MiFish.SWARM4.counts_LULU.txt
#worked!

#now for lulu
$R
> library(lulu)
> matchlist <- read.table("MiFish.SWARM4nc_seeds.obi.match_list.txt", header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)
> otutab <- read.csv("MiFish.SWARM4.counts_LULU.txt", sep='\t', header=TRUE, as.is=TRUE, row.names = 1)
> curated_result <- lulu(otutab, matchlist)
> write.csv(curated_result$curated_table,"MiFish.swarm4.LULU.curated.csv")

#############
#Now we count with two main table files: 1) a file containing the read numbers per sample for each MOTU generated during SWARM, COI_reefworms_swarm13_output.counts.csv, or the one corrected by LULU, COI_reefworms.LULU.curated.csv, and 2) a file containing the taxonomic information of the centroid sequences of each of the MOTUs, COI_reefworms.ecotag.annotated.csv (or COI_reefworms_taxonomy.csv if you used the BOLD querying strategy in Gate 4). Combining these two will generate a table with both the abundance in each sample and the taxonomic assignation of the MOTUs. But that table is not perfect and final, as there are some more refinements that will improve the presentation and analysis. At the beginning we will start with the line that included LULU curation, and the line that did NOT include it starts at 5.2 (black).
#5.1 First of all, since the counts and the taxonomy tables come from two different software, the formatting is different. We have to make them compatible, and it is just a matter of a little text editing.
sed -i -e 's/""/"id"/g' MiFish.swarm4.LULU.curated.csv 
sed -i -e 's/"//g' MiFish.swarm4.LULU.curated.csv
sed -i -e 's/,/;/g' MiFish.swarm4.LULU.curated.csv
#5.2 Now we combine both files with another one of the OwiTools scripts.
../owitools/owi_combine -i NCSharks_MiFish_results.csv -a MiFish.swarm4.LULU.curated.csv -o MiFish.results.combMOTUs.csv

ahhh_ok_I_think_the error is because the obi3_output_csv_has_;_at_the_end_of_each_line
sed 's/;//g' NCSharks_MiFish_results.csv > NCSharks_Mifish_results_for_combining.csv
../owitools/owi_combine -i NCSharks_Mifish_results_for_combining.csv -a MiFish.swarm4.LULU.curated.csv -o MiFish.results.combMOTUs.csv

the non-lulu_count_file is: this is used in the physalia course in this step.
MiFish.SWARM4nc_output.counts.csv

../owitools/owi_combine -i NCSharks_Mifish_results_for_combining.csv -a MiFish.SWARM4nc_output.counts.csv -o MiFish.results.combMOTUs.csv
stillll not working.
Reading ecotag database...
Ecotag database read including 104 total MOTUs.
Reading abundance database...
Abundances database read including 282 total MOTUs and 157 samples.
Error in fix.by(by.x, x) : 'by' must specify a uniquely valid column
Calls: merge -> merge.data.frame -> fix.by
Execution halted


or maybe I need to put the lulu results through obi3.

less MiFish.SWARM4.counts_LULU.txt
(eDNA) Eldridge-Wiselys-2022-MacbookPro:results Eldridge$ sed -i -e 's/id/""/g' MiFish.swarm4.LULU.curated.csv
(eDNA) Eldridge-Wiselys-2022-MacbookPro:results Eldridge$ less MiFish.swarm4.LULU.curated.csv
(eDNA) Eldridge-Wiselys-2022-MacbookPro:results Eldridge$ sed -i -e 's/""//g' MiFish.swarm4.LULU.curated.csv
(eDNA) Eldridge-Wiselys-2022-MacbookPro:results Eldridge$ less MiFish.swarm4.LULU.curated.csv
(eDNA) Eldridge-Wiselys-2022-MacbookPro:results Eldridge$ cat headers_MiFish.SWARM4.counts_LULU.txt MiFish.swarm4.LULU.curated.csv > MiFish.swarm4.LULU.curated_w_headers.csv

added headers to the lulu curated csv.

still didn't work

ok, here's the owi_combine from the physalia course
owi_combine -i ULOY.ecotag.fasta.annotated.csv -a ULOY.SWARM13nc_output.counts.csv -o ULOY.All_MOTUs.csv

need to check the format of the ecotag_db file: ULOY.ecotag.fasta.annotated.csv
###nope none of this below worked#####
obi 3 command instead of obi combine from obi2 would be:
#obi2:# owi_combine -i NCSharks_MiFish_results.csv -a MiFish.SWARM4nc_output.counts.csv -o FISH.All_MOTUs.csv
scp MiFish.swarm4.LULU.curated.csv eldridge@leopardus.snrenet.arizona.edu:/storage/eldridge/Database_for_fish_and_crustaceans/obitools/
scp ../02_obiconvert_to_fasta/MiFish.SWARM4nc_output.counts.csv eldridge@leopardus.snrenet.arizona.edu:/storage/eldridge/Database_for_fish_and_crustaceans/obitools/

obi export --tab-output fish/MiFish_swarm4_nosingles_tax_assigned > NCSharks_MiFish_results.fasta.csv
obi import --tabular-input MiFish.swarm4.LULU.curated.csv fish/LULU_result #failed and can't delete or overwrite it so,,
obi import --sep ';' --tabular-input MiFish.swarm4.LULU.curated.csv fish/LULU_result1
obi cat -c fish/MiFish_swarm4_nosingles_tax_assigned -c fish/LULU_result1 fish/combined_LULU_obi3
obi export --tab-output fish/combined_LULU_obi3 > MiFish_LULU_obi3.combined.tsv
#########################################


Maybe these files are what's needed as input for taxontabletools and I might not even need to owi_combine.

Nope, they would need lots of editing to get them into the right format.

First, I'm going to try vsearch with the CRABS database I made.
CRABS_db: FINAL_MiFish_taxonomy.sintax.fasta
query: MiFish.4seeds_nonsingleton.fasta

vsearch --sintax MiFish.4seeds_nonsingleton.fasta --db ../databases/FINAL_MiFish_taxonomy.sintax.fasta --tabbedout NC_Sharks_MiFish_sintax_taxo_assigned.csv --strand both
Classified 217 of 282 sequences (76.95%)
identified every query sequence to every taxonomic level!  but each level is followed by a numberinparentheses. Ahh, that number is the bootstrap confidence interval.

vsearch --sintax MiFish.4seeds_nonsingleton.fasta --db ../databases/FINAL_MiFish_taxonomy.sintax.fasta --tabbedout NC_Sharks_MiFish_sintax_taxo_assigned_97.csv --sintax_cutoff .97 --strand both
Classified 217 of 282 sequences (76.95%)


definitely have to LULU this.
I re-did LULU on the swarm 4 output files 
library(lulu)
matchlist <- read.table("MiFish.SWARM4nc_seeds.obi.match_list.txt", header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)
otutab <- read.csv("MiFish.SWARM4.counts_LULU.txt", sep='\t', header=TRUE, as.is=TRUE, row.names = 1)
curated_result <- lulu(otutab, matchlist)
write.csv(curated_result$curated_table,"MiFish.swarm4.LULU.curated.csv")

so this can be my otu table and the result of vsearch can be my taxontable for taxontabletools... with some formatting
but I just don't like sintax for taxonomic assignment as much as ecotag. Maybe the LULU curation helped (if I can figure out how to incorporate it) 
or an LCA script like that in eDNAFlow will help curate it more

I wish I could get ecotag to work outside of obitools3's weird and difficult to use system.

I guess editing the terrible formatting from obi3 is the way I want to go.

OK I guess we're at a crossroads. I could either mess with the formatting of the syntax taxonomy file or I could use the Obitools3 output into MetabaR. I saw where OB three has the output formats for metabaR so it seems like that's the way the pipeline wants me to go.
I really like the ecotags algorithm much better than the syntax algorithm for assigning taxonomy. I trust it more. So the first thing I will try is putting the obitools3 format into Metabar in R.

obi export --metabaR-output fish/MiFish_swarm4_nosingles_tax_assigned -o NCSharks_MiFish_results.metabar.tab
#Exception: Prefix needed when exporting for metabaR (--metabaR-prefix option)
obi export --metabaR-output fish/MiFish_swarm4_nosingles_tax_assigned --metabaR-prefix MiFish -o NCSharks_MiFish_results.metabar.tab
#KeyError: 'Cannot access to column sample in view MiFish_swarm4_nosingles_tax_assigned'


#There are no sample names in there so I will have to import a metadata file which is fine, if I knew what FORMAT to make it in!

I will try to reformat for TaxonTableTools before trying that.
In the meantime I will re-download the database to OBi3
#this will download all vertebrates except mammals:
wget -nH --cut-dirs=6 -A 'STD_VRT*.dat.gz' -R 'STD_MAM*.dat.gz' -m -np ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/std/
Wow! those 4 were really all!

downloading and installing metabaR
made a new folder in Savannah/NC_Sharks_obitools for the working directory.

formatted the desired 4 tables in Excel.  I used MiFish.swarm4.LULU.curated.csv, and NCSharks_MiFish_results.fasta (from obi3 tax_assigned) and the sample sheet for sequencing and the spreadsheets of the sample info(both sampling and lab)

I created the table PCRs by looking at the soil_euk example data and filling in the information from the sample sheet for sequencing
I created a new column in pcrs for sample_id and made it Shark# with the number corresponding to the NC fecal swab number when available and made new Shark#s for unmatched stomachs.  I put these in the Stomach Content Data table that Savannah sent me and put it in the shared google drive
I created another new column in pcrs for material and named them either fecal, stomach, or negative
I added columns for type sample or control and control_type the controls in this experiment were all PCR negatives
I changed IJKMNP of plate 2 (which were the names of the adapters, to ABCDEF)

I sorted Motus and Reads by Motu# order then copied and pasted reads (number of reads of each motu in each sample from LULU curated swarm4 file)
into a new excel sheet but used paste-special, transpose to get it so that the samples were columns and motus were rows.
no first cell label on either of these is needed (don't use one)

Importing caused some issues if the files weren't copied and pasted into new sheets because of empty row or column names/cells at the ends of the data.

Eventually, it worked!

Now, continuing with the tutorial! https://metabarfactory.github.io/metabaR/articles/metabaRF-vignette.html#tutorial-with-the-soil_euk-dataset

12/5/22 got to the end of the metabar analysis pipeline with the MiFish reads, saved all the input and output files from Metabar in /Users/Eldridge/Desktop/Savannah/NC_Sharks_obitools/2nd_try_demuxed_cleaned/MetabaR

OK.  Now for the crustaceans.  This document has gotten long and confusing, so for the sake of clarity,I'm going to make a new document for the crustaceans, doing exactly the same steps that worked here for the fish.
#haha, that would be too easy.  Because the crustaceans were giving me trouble, I brought them back over here so that FINAL_NC_Sharks_bioinformatics... remains clean.

#NC_Sharks_Crustaceans-only

while read line; do    
AdapterRemoval --threads 10 --file1 "../BerryCrust/""$line""_R1_BerryCrust.fastq" --file2 "../BerryCrust/""$line""_R2_BerryCrust.fastq" \
                    --collapse --trimns --trimqualities \
                    --minquality 20 \
                    --minalignmentlength 12 \
                    --basename "$line""_BerryCrust"; done < ../Sample_list.txt

mkdir ../01b_fastqc_merged_filtered
fastqc *.collapsed -o ../01b_fastqc_merged_filtered/
cd ../01b_fastqc_merged_filtered/
multiqc .

# I think they all need adapter removal done because the multiqc still says many of them have Nextera transposase adapter contamination.  The lengths look really good though.
# I checked some of the files of adapter contamination concern from fastqc/multiqc and both the MiFish and BerryCrust reads still have the forward PCR1 primer adapter TCGTCGGCAGCGTC and often followed by the 3' end of the illumina adapter AGATGTGTATAAGAGACAG hanging out on the ends of the reads.
      
# I used cutadapt again to clean the reads using the fwd and revcomp reads of each section of adapter applied in PCR1 and 2 other than the MiFish and BerryCrust primers      

while read line; do    
cutadapt --times 10 -b TCGTCGGCAGCGTC -b AGATGTGTATAAGAGACAG -b GTCTCGTGGGCTCGG -b AATGATACGGCGACCACCGAGATCTACAC -b CAAGCAGAAGACGGCATACGAGAT -b GACGCTGCCGACGA -b CTGTCTCTTATACACATCT -b CCGAGCCCACGAGAC -b GTGTAGATCTCGGTGGTCGCCGTATCATT -b ATCTCGTATGCCGTCTTCTGCTTG -o "../01c_cutadapt_merged_qualityfiltered/""$line""_BerryCrust.collapsed.AR" "$line""_BerryCrust.collapsed" ; done < ../Sample_list.txt
#re-check for adapter contamination

cd ../01c_cutadapt_merged_qualityfiltered/
fastqc * -o ../01d_fastqc_merged_filtered_adapter_removed/
cd ../01d_fastqc_merged_filtered_adapter_removed
multiqc .

#Looks much better!


#Add Annotations

for file in *_BerryCrust.collapsed.AR; do
sample=${file#*dir/}
sample=${sample%_BerryCrust.collapsed.AR}
obiannotate -S sample:${sample} -S assay:BerryCrust --length ${sample}_BerryCrust.collapsed.AR > ${sample}_BerryCrust.collapsed.AR.annotated.fastq
done

#Now I'd like to filter for sequences between 140-260 bp long.  Also keep quality over 20, maximum 0 Ns and 1 expected errors.
mkdir ../01e_trimmed_QF
mkdir ../01e_trimmed_QF/discards
while read line; do
cutadapt -q 20 --max-n 0 --max-ee 1 -m 140 -M 260 --too-short-output "../01e_trimmed_QF/discards/""$line""_BerryCrust.under140bp" --too-long-output "../01e_trimmed_QF/discards/""$line""_BerryCrust.over260bp" -o "../01e_trimmed_QF/""$line""_BerryCrust.final.annotated.fastq" "$line""_BerryCrust.collapsed.AR.annotated.fastq"; done < ../Sample_list.txt
#check for length, quality, and no contamination
cd ../01e_trimmed_QF/
fastqc * -o ../01f_fastqc_trimmed_QF_final_annotated/
cd ../01f_fastqc_trimmed_QF_final_annotated
multiqc .

#Looks great!  they're ready for the next steps now

cd /Users/Eldridge/Desktop/Savannah/NC_Sharks_obitools/2nd_try_demuxed_cleaned
scp -r eldridge@leopardus.snrenet.arizona.edu:/storage/eldridge/NC_Sharks_data/demuxed_primers_w_cutadapt/following_eDNAFlow_w_demuxed_data/01e_trimmed_QF/ .
mkdir 02_obiconvert_to_fasta

while read line; do  
obiconvert --sanger --fasta-output \
"$line""_BerryCrust.final.annotated.fastq" > "../02_obiconvert_to_fasta/""$line""_BerryCrust""".fasta""; done < ../Sample_list.txt

#11/14/2022 in obitools environment on my new computer in 02_obiconvert_to_fasta directory
#concatenate into one file per primer set and count sequences per sample.
cat *_BerryCrust.fasta  > BerryCrust.allsamples.fasta
obistat -c sample -a seq_length BerryCrust.allsamples.fasta > sample_stats_BerryCrust.txt 

**********************
obisplit -t sample BerryCrust.allsamples.fasta 
while read line; do  
rm "$line"".fasta"; done < ../Sample_list.txt

#only keep unique sequences in each sample
for i in *_BerryCrust.fasta ; do obiuniq $i >  ${i/.fasta/.unique.fasta}; done

#change from obitools format to vsearch format using owitools scripts
cd ..
mkdir owitools
cd owitools
cp -r ~/Desktop/Metabarcoding_Physalia/ubuntu/owi_tools/ .
chmod 777 *
which RScript 
result on my computer: /usr/local/bin/RScript
copy and paste that into the first line of the owi_script that I want to use
also had to go into RStudio and install.packages("optparse")

for i in *_BerryCrust.unique.fasta ; do ./../owitools/owi_obisample2vsearch -i $i; done
for i in *_BerryCrust.unique.vsearch.fasta ; do vsearch --uchime_denovo $i --sizeout --minh 0.90 --nonchimeras ${i/.unique.vsearch.fasta/.nonchimeras.fasta} --chimeras ${i/.unique.vsearch.fasta/.chimeras.fasta} --uchimeout ${i/.unique.vsearch.fasta/.uchimeout.txt}; done


#Concatenating non chimaeras in a single file
*_BerryCrust.nonchimeras.fasta > BerryCrust.nonchimeras.fasta

#Returning to obitools format
../owitools/owi_vsearch2obifasta -i BerryCrust.nonchimeras.fasta
#Dereplicating sequences in MiFish.nonchimeras.fasta file

obiuniq -m sample BerryCrust.nonchimeras.vsearch.fasta > BerryCrust.unique.fasta
#change names of motus to a short more informative identifier

obiannotate --seq-rank BerryCrust.unique.fasta | obiannotate --set-identifier '"BerryCrust_%07d" % seq_rank' > BerryCrust.new.fasta
#obtain a table of abundances

obitab -o BerryCrust.new.fasta >  BerryCrust.new.tab

#back to vsearch format

#../owitools/owi_obifasta2vsearch -i BerryCrust.new.fasta -o BerryCrust.vsearch.fasta #not working right

# I copied a new version of owi_obifasta2vsearch from the original folder and modified it to know that the name was the first 18 characters instead of 14 characters for BerryCrust_%07d




#cluster with swarm :BerryCrust
conda activate eDNA
swarm -d 1 -f -z -t 40 -o BerryCrust.SWARM1nc_output -s BerryCrust.SWARM1nc_stats -w BerryCrust.SWARM1nc_seeds.fasta BerryCrust.vsearch.fasta

Duplicated sequence identifier: BerryCrust_000
#hmm now it seems to want only 3 digits.  I hope (but am pretty sure) that's enough
less BerryCrust.new.fasta
tail BerryCrust.new.fasta #156805 unique sequences
needs at least 6 digits, I'll do 7

obiannotate --seq-rank BerryCrust.unique.fasta | obiannotate --set-identifier '"BerryCrust_%07d" % seq_rank' > BerryCrust.new.fasta
# I copied a new version of owi_obifasta2vsearch from the original folder and modified it to know that the name was the first 18 characters instead of 14 characters for BerryCrust_%07d


swarm -d 2 -z -t 40 -o BerryCrust.SWARM2nc_output -s BerryCrust.SWARM2nc_stats -w BerryCrust.SWARM2nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 3 -z -t 40 -o BerryCrust.SWARM3nc_output -s BerryCrust.SWARM3nc_stats -w BerryCrust.SWARM3nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 4 -z -t 40 -o BerryCrust.SWARM4nc_output -s BerryCrust.SWARM4nc_stats -w BerryCrust.SWARM4nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 5 -z -t 40 -o BerryCrust.SWARM5nc_output -s BerryCrust.SWARM5nc_stats -w BerryCrust.SWARM5nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 6 -z -t 40 -o BerryCrust.SWARM6nc_output -s BerryCrust.SWARM6nc_stats -w BerryCrust.SWARM6nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 7 -z -t 40 -o BerryCrust.SWARM7nc_output -s BerryCrust.SWARM7nc_stats -w BerryCrust.SWARM7nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 8 -z -t 40 -o BerryCrust.SWARM8nc_output -s BerryCrust.SWARM8nc_stats -w BerryCrust.SWARM8nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 9 -z -t 40 -o BerryCrust.SWARM9nc_output -s BerryCrust.SWARM9nc_stats -w BerryCrust.SWARM9nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 10 -z -t 40 -o BerryCrust.SWARM10nc_output -s BerryCrust.SWARM10nc_stats -w BerryCrust.SWARM10nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 11 -z -t 40 -o BerryCrust.SWARM11nc_output -s BerryCrust.SWARM11nc_stats -w BerryCrust.SWARM11nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 12 -z -t 40 -o BerryCrust.SWARM12nc_output -s BerryCrust.SWARM12nc_stats -w BerryCrust.SWARM12nc_seeds.fasta BerryCrust.vsearch.fasta



#R-script
#d values for swarm test y=number of swarms MiFish
# x=d-values y=number of swarms
x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
y<- c(5712,3924,1676,977,713,572,485,436,386,360,336,304,260,237)
plot(x,y)
xtrim <-c(3,4,5,6,7,8,9,10,11,12,13,14)
ytrim<-c(1676,977,713,572,485,436,386,360,336,304,260,237)
plot(xtrim,ytrim)

#d values for swarm test y=number of swarms BerryCrust
a <- c(1,2,3,4,5,6,7,8,9,10,11,12)
b<- c(5937,4687,2458,1699,1337,1151,1063,1008,966,925,898,862)
plot(a,b)
atrim <-c(3,4,5,6,7,8,9,10,11,12)
btrim<-c(2458,1699,1337,1151,1063,1008,966,925,898,862)
plot(atrim,btrim)

#BerryCrust results
Kept only 656 clusters of size greater than or equal to 2 reads.
Reading tabulated database. This could take a while...
Database read including 156805 total different sequences and 157 samples.
Kept only 0 sequences for calculations.         #this is a problem.  It's the same for d=1-12
File BerryCrust.SWARM5nc_output.counts.csv written

#Remove singletons

 ../owitools/owi_vsearch2obifasta -i BerryCrust.SWARM5nc_seeds.fasta -o BerryCrust.SWARM5nc_seeds.obi.fasta
obigrep -p 'count>1' BerryCrust.SWARM5nc_seeds.fasta > BerryCrust.seeds_nonsingleton.fasta
#no nonsingletons left at d=5 swarm

../owitools/owi_recount_swarm BerryCrust.SWARM4nc_output BerryCrust.swarm4.new.tab
# did the same for 1-12, none had any non-singletons left

#so I processed the reads (uniques,non-chimaeras) in dada2 on my computer and exported them to a new fasta file
#then in a new folder Savannah/dada2_to_obi3
obiannotate --seq-rank ASVs_without_bimeras_BerryCrust.fasta | obiannotate --set-identifier '"BerryCrust_%07d" % seq_rank' > BerryCrust.new.fasta
obitab -o BerryCrust.new.fasta >  BerryCrust.new.tab


#took results of dada2 following french Abyss tutorial and tried swarm on them.
on my computer:
conda activate eDNA

conda update swarm
swarm version=3.13
swarm -d 1 -f -z -t 40 -o BerryCrust.dada_SWARM1nc_output -s BerryCrust.dada_SWARM1nc_stats -w BerryCrust.dada_SWARM1nc_seeds.fasta 01e_trimmed_QFBerryCrust_ASVs.fasta
#says it needs the count information in the headers which can be gotten by dereplicating the sequences.  



#obi3 for taxonomic assignment#
#on leopardus in obitools3 

conda activate obi3
source obi3-env/bin/activate

#for BerryCrust, I only downloaded invertebrates (I left in the commands to not include humans or environmental sequences)
#in /storage/eldridge/Database_for_fish_and_crustaceans/obitools/crustaceans_EMBL$ 
wget -nH --cut-dirs=6 -A 'STD_INV*.dat.gz' -R 'STD_HUM*.dat.gz','STD_ENV*.dat.gz' -m -np ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/std/
#this downloaded 9 files

obi import --embl EMBL crust/Berry_refs

wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

obi import --taxdump taxdump.tar.gz crust/taxonomy/my_tax



#######BECAUSE THE CRUSTACEANS WON'T PLAY NICE WITH SWARM, WE'RE GOING BACK TO DADA2 USING THE FRENCH ABYSS PROJECT SCRIPTS#######

first make an rdp formatted taxonomy database in CRABS.
Savannah did everything except the new pga step, so I'll start from there:
in /Users/Eldridge/Desktop/Savannah/CRABS/BerryCrust_DB
crabs pga --input BerryCrust_db_merged.fasta --output BerryCrust_db_pga.fasta --database BerryCrust_db_PCR.fasta --fwd GGGACGATAAGACCCTATA  --rev ATTACGCTGTTATCCCTAAAG  --speed medium --percid 0.95 --coverage 0.95 --filter_method strict
crabs assign_tax -i BerryCrust_db_pga.fasta -o BerryCrust_db_tax.tsv -a nucl_gb.accession2taxid -t nodes.dmp -n names.dmp -w “yes” 
crabs dereplicate --input BerryCrust_db_tax.tsv --output BerryCrust_db_dereplicate.tsv --method uniq_species
crabs seq_cleanup --input BerryCrust_db_dereplicate.tsv --output BerryCrust_db_cleanup.tsv --minlen 70 --maxlen 270 --maxns 0 --species yes --nans 1
crabs tax_format --input BerryCrust_db_cleanup.tsv --output FINAL_BerryCrust_taxonomy.fasta --format rdp

#I'm Not liking that insects are included in the database.  I'm going to re-make it
crabs db_download --source ncbi --database nucleotide --query ‘"Crustacea"[Organism] OR "Decapoda"[Organism] AND 16S[All Fields] AND animals[filter] AND biomol_genomic[PROP] AND ddbj_embl_genbank[filter] AND is_nuccore[filter] AND mitochondrion[filter])’ --output Crustaceans_and_Decapods_16S_ncbi.fasta --keep_original yes --email eldridgewisely@arizona.edu --batchsize 5000

#actually, it's fine, perhaps I can just remove insects later in the script because at least one of the shrimp iD'ed in dada2 rdp taxoassignment is not listed as decapod or crustacean
********************

#ok, starting with the french abyss project r-script dada2main.R downloaded from: and stored in Savannah/following_dada2_method
cd /Users/Eldridge/Desktop/Savannah/following_dada2_method/02_cleaned_trimmed/minlen2_dec11-22

#these are the commands that were used in November for these reads except with --minimum-length 2 added
conda activate eDNA

while read line; do    
cutadapt --revcomp --times 10 -b TCGTCGGCAGCGTC -b AGATGTGTATAAGAGACAG -b GTCTCGTGGGCTCGG -b AATGATACGGCGACCACCGAGATCTACAC -b CAAGCAGAAGACGGCATACGAGAT --minimum-length 2 -o "$line""_BerryCrust_R1.cleaned.fastq" "../../NC_Sharks_cutadapt/""$line""_R1_BerryCrust.fastq"; done < ../../Sample_list.txt
while read line; do
cutadapt --revcomp --times 10 -b GACGCTGCCGACGA -b CTGTCTCTTATACACATCT -b CCGAGCCCACGAGAC -b GTGTAGATCTCGGTGGTCGCCGTATCATT -b ATCTCGTATGCCGTCTTCTGCTTG --minimum-length 2 -o "$line""_BerryCrust_R2.cleaned.fastq" "../../NC_Sharks_cutadapt/""$line""_R2_BerryCrust.fastq"; done < ../../Sample_list.txt

while read line; do    
cutadapt --revcomp --times 10 -b TCGTCGGCAGCGTC -b AGATGTGTATAAGAGACAG -b GTCTCGTGGGCTCGG -b AATGATACGGCGACCACCGAGATCTACAC -b CAAGCAGAAGACGGCATACGAGAT --minimum-length 2 -o "$line""_MiFish_R1.cleaned.fastq" "../../NC_Sharks_cutadapt/""$line""_R1_MiFish.fastq"; done < ../../Sample_list.txt
while read line; do
cutadapt --revcomp --times 10 -b GACGCTGCCGACGA -b CTGTCTCTTATACACATCT -b CCGAGCCCACGAGAC -b GTGTAGATCTCGGTGGTCGCCGTATCATT -b ATCTCGTATGCCGTCTTCTGCTTG --minimum-length 2 -o "$line""_MiFish_R2.cleaned.fastq" "../../NC_Sharks_cutadapt/""$line""_R2_MiFish.fastq"; done < ../../Sample_list.txt

cutadaptversion 4.1 python 3.8.2

#checking them: 

conda activate qc
fastqc *.fastq -o ../../03_qc_check/minlen2_dec11-22/
in R, I'm trying to process both at the same time so I haveto usea different separator betweenb the sampleandassay.
rename -s _R .R *.fastq


#using merged reads like the Fish
CRABS Eldridge$ crabs tax_format --input BerryCrust_db_cleanup.tsv --output FINAL_BerryCrust_taxonomy.dad.fasta --format dad
same for dads and rdp
in /Users/Eldridge/Desktop/Savannah/CRABS

#sent preliminary results (merged and cleaned like MiFish, then put into dada2 for ASVs instead of swarm for OTUs, then rdp taxonomic assignment in dada2 instead of obi3) to Savannah

#re-doing the demultiplexing of NC7 in /Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/NC_Sharks_cutadapt
cutadapt -e 0.15 --no-indels -g MiFish=GTCGGTAAAACTCGTGCCAGC -g BerryCrust=GGGACGATAAGACCCTATA -G MiFishUR=CATAGTGGGGTATCTAATCCCAGTTTG -G BerryCrustR=ATTACGCTGTTATCCCTAAAG -o NC7_R1_{name}.fastq -p NC7_R2_{name}.fastq ../../../ORIGINAL\ NC_Sharks_data/NC7xx-B1_CTAGGTGA-GCATAACG_S79_L001_R1_001.fastq.gz ../../../ORIGINAL\ NC_Sharks_data/NC7xx-B1_CTAGGTGA-GCATAACG_S79_L001_R2_001.fastq.gz
fixed the empty NC7 problem!

#next,I will re-do the adapter removal with adapterremoval
#AdapterRemoval --file1 reads_1.fq --file2 reads_2.fq --basename output_paired --trimns --trimqualities --collapse
#nope not yet, I think I see the problem with the way I did cutadapt.  I didn't treat them as paired-end reads when removing the adapters.

cd /Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/clean_w_cutadapt_dec12-22
while read line; do
cutadapt -q 20 --trim-n --minimum-length 10 -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CCGAGCCCACGAGAC -a ATCTCGTATGCCGTCTTCTGCTTG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A GACGCTGCCGACGA -A GTGTAGATCTCGGTGGTCGCCGTATCATT -e .1 -o "$line""_MiFish_R1.cleaned.fastq" -p "$line""_MiFish_R2.cleaned.fastq" "../NC_Sharks_cutadapt/""$line""_R1_MiFish.fastq" "../NC_Sharks_cutadapt/""$line""_R2_MiFish.fastq"; done < ../Sample_list.txt

mkdir fastqc
fastqc * -o fastqc/
cd fastqc
multiqc .

fastqc says it's finding Nextera transposase sequences in the reverse reads.

I'll try trimmomatic
conda install trimmomatic (in qc)
in /Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/clean_w_trimmomatic_dec12-22
while read line; do
trimmomatic PE -trimlog trimmomatic.log "../NC_Sharks_cutadapt/""$line""_R1_MiFish.fastq" "../NC_Sharks_cutadapt/""$line""_R2_MiFish.fastq" "$line""_MiFish_R1.trimmed.fastq"  "$line""_MiFish_R1.discards.fastq" "$line""_MiFish_R2.trimmed.fastq" "$line""_MiFish_R2.discards.fastq" -validatePairs ILLUMINACLIP:../NexteraPE-PE.fa:2:30:10:2:true SLIDINGWINDOW:4:20 MINLEN:10;done < ../Sample_list.txt

mkdir fastqc
fastqc * -o fastqc/
cd fastqc
multiqc .

There still appears to be adapter contamination above 200bp in half of the reads
while read line; do
trimmomatic PE -trimlog trimmomatic.log "../NC_Sharks_cutadapt/""$line""_R1_MiFish.fastq" "../NC_Sharks_cutadapt/""$line""_R2_MiFish.fastq" "$line""_MiFish_R1.trimmed.fastq"  "$line""_MiFish_R1.discards.fastq" "$line""_MiFish_R2.trimmed.fastq" "$line""_MiFish_R2.discards.fastq" -validatePairs ILLUMINACLIP:../NexteraPE-PE.fa:2:30:10:2:true SLIDINGWINDOW:4:15 CROP:200 MINLEN:30;done < ../Sample_list.txt

#edited the NexteraPE-PE.fa file to do palindrome (read-through) trimming of everything.  saved as NexteraPE-PE2.fa
while read line; do trimmomatic PE -trimlog trimmomatic.log "../NC_Sharks_cutadapt/""$line""_R1_MiFish.fastq" "../NC_Sharks_cutadapt/""$line""_R2_MiFish.fastq" "$line""_MiFish_R1.trimmed.fastq"  "$line""_MiFish_R1.discards.fastq" "$line""_MiFish_R2.trimmed.fastq" "$line""_MiFish_R2.discards.fastq" -validatePairs ILLUMINACLIP:../NexteraPE-PE2.fa:2:30:10:2:true SLIDINGWINDOW:4:15 CROP:225 MINLEN:30;done < ../Sample_list.txt

still issues with this.  I will try combining the two.  
mkdir trimmed (in /Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/clean_w_cutadapt_dec12-22)
cd trimmed

while read line; do
trimmomatic PE -trimlog trimmomatic.log "../""$line""_MiFish_R1.cleaned.fastq" "../""$line""_MiFish_R2.cleaned.fastq" "$line""_MiFish_R1.trimmed.fastq"  "$line""_MiFish_R1.discards.fastq" "$line""_MiFish_R2.trimmed.fastq" "$line""_MiFish_R2.discards.fastq" -validatePairs ILLUMINACLIP:../../NexteraPE-PE.fa:2:30:10:2:true SLIDINGWINDOW:4:15;done < ../../Sample_list.txt

#trimmomatic caught a few and I like the sliding window trimming but there's still supposedly adapter contamination in some of them.
so I added the only other sequences that were part of my primers to the NexteraPE file (even though I included the rc of them in the cutadapt command with where they were supposed to be if they showed up.)

while read line; do
trimmomatic PE -trimlog trimmomatic.log "../""$line""_MiFish_R1.cleaned.fastq" "../""$line""_MiFish_R2.cleaned.fastq" "$line""_MiFish_R1.trimmed.fastq"  "$line""_MiFish_R1.discards.fastq" "$line""_MiFish_R2.trimmed.fastq" "$line""_MiFish_R2.discards.fastq" -validatePairs ILLUMINACLIP:../../NexteraPE-PE2.fa:2:30:10:2:true SLIDINGWINDOW:4:15;done < ../../Sample_list.txt

when these sequences show up in the cleaned files, it's at the beginning... 5th character etc.  mostly in the r2s... so a right trim of r2 might be in order(-G in cutadapt)


Dec-15-22
#cutadapt cleaning, adding MiFish UR and UF and rc of each and the 3' end of the illumina adapter added to both, not just the rc of it on the 3' end like before
cd Desktop/Savannah/NC_Sharks_data/Cutadapt
mkdir clean_w_cutadapt_dec15-22
cd clean_w_cutadapt_dec15-22

while read line; do
cutadapt -q 20 --trim-n --minimum-length 10 -g AGATGTGTATAAGAGACAG -a CAAACTGGGATTAGATACCCCACTATG -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CCGAGCCCACGAGAC -a ATCTCGTATGCCGTCTTCTGCTTG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A GACGCTGCCGACGA -A GTGTAGATCTCGGTGGTCGCCGTATCATT -A GCTGGCACGAGTTTTACCGACNNNN -G AGATGTGTATAAGAGACAG -G CATAGTGGGGTATCTAATCCCAGTTTG -e .1 -o "$line""_MiFish_R1.cleaned.fastq" -p "$line""_MiFish_R2.cleaned.fastq" "../NC_Sharks_cutadapt/""$line""_R1_MiFish.fastq" "../NC_Sharks_cutadapt/""$line""_R2_MiFish.fastq"; done < ../Sample_list.txt


these were added to the previous command:

AGATGTGTATAAGAGACAG 3' end of Illumina adapter (nextera) found in R reads on the left and some F reads on the left too. -g -G
GCTGGCACGAGTTTTACCGACNNNN rc of MiFishUF found in MiFish R reads on the right -A
CATAGTGGGGTATCTAATCCCAGTTTG MiFishUR found in left side of R reads usually attached to rc of MiFishUF -G
CAAACTGGGATTAGATACCCCACTATG rc of MiFish UR on the right side of F reads -a

GTCGGTAAAACTCGTGCCAGC MiFishUF-not seen in any reads (probably because it was used in demultiplexing)

#ok, I'll do the same for the crustaceans: here's the first command with all of the adapter sequences they have in common mostly found rc in the opposite read

while read line; do
cutadapt -q 20 --trim-n --minimum-length 10 -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CCGAGCCCACGAGAC -a ATCTCGTATGCCGTCTTCTGCTTG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A GACGCTGCCGACGA -A GTGTAGATCTCGGTGGTCGCCGTATCATT -e .1 -o "$line""_MiFish_R1.cleaned.fastq" -p "$line""_MiFish_R2.cleaned.fastq" "../NC_Sharks_cutadapt/""$line""_R1_MiFish.fastq" "../NC_Sharks_cutadapt/""$line""_R2_MiFish.fastq"; done < ../Sample_list.txt

adding 3' end of Illumina adapter (nextera) found in R reads on the left and some F reads on the left too.
-g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG
BerryCrustR on left of R reads 
-G NNNATTACGCTGTTATCCCTAAAG
rc of BerryCrustR on the right side of R1s(f)
-a CTTTAGGGATAACAGCGTAATNNN
rc of BerryCrustF on the right side of the R2s(r)
-A TATAGGGTCTTATCGTCCCNNN

GGGACGATAAGACCCTATA BerryCrustF -not seen in any reads (probably because it was used in demultiplexing)


while read line; do
cutadapt -q 20 --trim-n --minimum-length 10 -a CTTTAGGGATAACAGCGTAATNNN -A TATAGGGTCTTATCGTCCCNNN -G NNNATTACGCTGTTATCCCTAAAG -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CCGAGCCCACGAGAC -a ATCTCGTATGCCGTCTTCTGCTTG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A GACGCTGCCGACGA -A GTGTAGATCTCGGTGGTCGCCGTATCATT -e .1 -o "$line""_BerryCrust_R1.cleaned.fastq" -p "$line""_BerryCrust_R2.cleaned.fastq" "../NC_Sharks_cutadapt/""$line""_R1_BerryCrust.fastq" "../NC_Sharks_cutadapt/""$line""_R2_BerryCrust.fastq"; done < ../Sample_list.txt

#BerryCrust looks pretty good!  I think that perhaps adding the primer-based cleaning at the beginning instead of the end might be helping.
mkdir trimmed
cd trimmed
while read line; do
trimmomatic PE -trimlog trimmomatic.log "../""$line""_BerryCrust_R1.cleaned.fastq" "../""$line""_BerryCrust_R2.cleaned.fastq" "$line""_BerryCrust_R1.trimmed.fastq"  "$line""_BerryCrust_R1.discards.fastq" "$line""_BerryCrust_R2.trimmed.fastq" "$line""_BerryCrust_R2.discards.fastq" -validatePairs ILLUMINACLIP:../../NexteraPE-PE.fa:2:30:10:2:true SLIDINGWINDOW:4:15;done < ../../Sample_list.txt


I'll do that for the MiFish now:
first command:

#while read line; do
#cutadapt -q 20 --trim-n --minimum-length 10 -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CCGAGCCCACGAGAC -a ATCTCGTATGCCGTCTTCTGCTTG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A GACGCTGCCGACGA -A GTGTAGATCTCGGTGGTCGCCGTATCATT -e .1 -o "$line""_MiFish_R1.cleaned.fastq" -p "$line""_MiFish_R2.cleaned.fastq" "../NC_Sharks_cutadapt/""$line""_R1_MiFish.fastq" "../NC_Sharks_cutadapt/""$line""_R2_MiFish.fastq"; done < ../Sample_list.txt

primer-based cleaning added to the front of the command
while read line; do
cutadapt -q 20 --trim-n --minimum-length 10 -a CAAACTGGGATTAGATACCCCACTATG -G NNNNCATAGTGGGGTATCTAATCCCAGTTTG -A GCTGGCACGAGTTTTACCGACNNNN -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CCGAGCCCACGAGAC -a ATCTCGTATGCCGTCTTCTGCTTG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A GACGCTGCCGACGA -A GTGTAGATCTCGGTGGTCGCCGTATCATT -e .1 -o "$line""_MiFish_R1.cleaned.fastq" -p "$line""_MiFish_R2.cleaned.fastq" "../NC_Sharks_cutadapt/""$line""_R1_MiFish.fastq" "../NC_Sharks_cutadapt/""$line""_R2_MiFish.fastq"; done < ../Sample_list.txt

while read line; do
trimmomatic PE -trimlog trimmomatic.log "../""$line""_MiFish_R1.cleaned.fastq" "../""$line""_MiFish_R2.cleaned.fastq" "$line""_MiFish_R1.trimmed.fastq"  "$line""_MiFish_R1.discards.fastq" "$line""_MiFish_R2.trimmed.fastq" "$line""_MiFish_R2.discards.fastq" -validatePairs ILLUMINACLIP:../../NexteraPE-PE.fa:2:30:10:2:true SLIDINGWINDOW:4:15;done < ../../Sample_list.txt


######the below did nothing for the tiny amount of remaining adapter detected so using the .trimmed.fastq files for future.  Deleting the results of the below to save space.
while read line; do
cutadapt -q 20 --trim-n --minimum-length 10 -b AGATGTGTATAAGAGACAG -B AGATGTGTATAAGAGACAG -e .1 -o "$line""_MiFish_R1.final_cleaned.fastq" -p "$line""_MiFish_R2.final_cleaned.fastq" "../""$line""_MiFish_R1.trimmed.fastq" "../""$line""_MiFish_R2.trimmed.fastq"; done < ../../../Sample_list.txt
#looks like the final_clean didn't do anything about the remaining adapters according to multiqc


#for berryCrust
while read line; do
cutadapt -q 20 --trim-n --minimum-length 10 -b AGATGTGTATAAGAGACAG -B AGATGTGTATAAGAGACAG -e .1 -o "$line""_BerryCrust_R1.final_cleaned.fastq" -p "$line""_BerryCrust_R2.final_cleaned.fastq" "../""$line""_BerryCrust_R1.trimmed.fastq" "../""$line""_BerryCrust_R2.trimmed.fastq"; done < ../../../Sample_list.txt


###############try to re-do BerryCrust the FINAL_NC_Sharks_bioinformatics_code.txt way########### still no sequences left for BerryCrust after swarm
I guess I'll collapse the BerryCrust trimmed reads for this analysis.(but I can also use them as-is in dada2 or obitools3)
#############
BerryCrust Obitools3 on Panthera

#first annotate and concatenate the cleaned reads
on my computer with the old obitools installed:
while read line; do
obiannotate -S sample:${line} -S assay:BerryCrust --length ${line}_BerryCrust_R1.trimmed.fastq > ${line}_BerryCrust_R1.trimmed.annotated.fastq; done < ../../Sample_list.txt

while read line; do
obiannotate -S sample:${line} -S assay:BerryCrust --length ${line}_BerryCrust_R2.trimmed.fastq > ${line}_BerryCrust_R2.trimmed.annotated.fastq; done < ../../Sample_list.txt


cat *_BerryCrust_R2.trimmed.annotated.fastq > NC-BerryCrust_annotated_R2.fastq
cat *_BerryCrust_R1.trimmed.annotated.fastq > NC-BerryCrust_annotated_R1.fastq

#wolf tutorial commands... adapted on panthera and noted on there.  These are just a guide.
obi import crust_tutorial/crust_F.fastq crust/reads1
obi import crust_tutorial/crust_R.fastq crust/reads2

obi alignpairedend -R crust/reads2 crust/reads1 crust/aligned_reads

obi stats -a score_norm crust/aligned_reads

obi grep -p "sequence['score_norm'] > 0.8" crust/aligned_reads crust/good_sequences

obi uniq -m sample crust/identified_sequences crust/dereplicated_sequences

obi annotate -k COUNT -k MERGED_sample crust/dereplicated_sequences crust/cleaned_metadata_sequences

obi grep -p "len(sequence)>=80 and sequence['COUNT']>=10" crust/cleaned_metadata_sequences crust/denoised_sequences

obi clean -s MERGED_sample -r 0.05 -H crust/denoised_sequences crust/cleaned_sequences

#obi import --embl EMBL crust/embl_refs 
#or


#obi import v05_refs.fasta.gz crust/v05_refs


#wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
#md5sum taxdump.tar.gz compare with https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5

obi import --taxdump taxdump.tar.gz crust/taxonomy/my_tax

obi ecopcr -e 3 -l 50 -L 150 -F TTAGATACCCCACTATGC -R TAGAACAGGCTCCTCTAG --taxonomy crust/taxonomy/my_tax crust/embl_refs crust/v05_refs

obi grep --require-rank=species --require-rank=genus --require-rank=family --taxonomy crust/taxonomy/my_tax crust/v05_refs crust/v05_refs_clean

obi build_ref_db -t 0.97 --taxonomy crust/taxonomy/my_tax crust/v05_refs_uniq_clean crust/v05_db_97

obi ecotag -m 0.97 --taxonomy crust/taxonomy/my_tax -R crust/v05_db_97 crust/cleaned_sequences crust/assigned_sequences

obi stats -c SCIENTIFIC_NAME crust/assigned_sequences

obi align -t 0.95 crust/assigned_sequences crust/aligned_assigned_sequences

obi history crust

obi history -d crust > crust.dot
obi history -d crust/cleaned_sequences > crust_one_view.dot

dot -Tx11 crust.dot

dot -Tpng crust.dot -o crust.png
open crust.png &

obi export --fasta-output crust/assigned_sequences -o crust_results.fasta

obi export --tab-output crust/aligned_assigned_sequences > crust_results.csv

#### my way, above, to see if the cleaning made a difference#####

cd /Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/BC_clean_w_cutadapt_dec15-22/trimmed
cd ../..
mkdir BC_my_pipeline_dec-16-22
cd BC_my_pipeline_dec-16-22/
mkdir merged
cd merged/
while read line; do    
AdapterRemoval --threads 10 --file1 "../../BC_clean_w_cutadapt_dec15-22/trimmed/""$line""_BerryCrust_R1.trimmed.fastq" --file2 "../../BC_clean_w_cutadapt_dec15-22/trimmed/""$line""_BerryCrust_R2.trimmed.fastq" \
                    --collapse --trimns --trimqualities \
                    --minquality 20 \
                    --minalignmentlength 12 \
                    --basename "$line""_BerryCrust"; done < ../../Sample_list.txt

mkdir ../fastqc_merged_filtered
fastqc *.collapsed -o ../fastqc_merged_filtered/
cd ../fastqc_merged_filtered/
multiqc .

#looks great!

for file in *_BerryCrust.collapsed; do
sample=${file#*dir/}
sample=${sample%_BerryCrust.collapsed}
obiannotate -S sample:${sample} -S assay:BerryCrust --length ${sample}_BerryCrust.collapsed > ${sample}_BerryCrust.collapsed.annotated.fastq
done

#filter for sequences between 100-260 bp long.  Also keep quality over 20, maximum 0 Ns and 1 expected errors.
mkdir merged_trimmed_QF
mkdir merged_trimmed_QF/discards
while read line; do
cutadapt -q 20 --max-n 0 --max-ee 1 -m 100 -M 260 --too-short-output "merged_trimmed_QF/discards/""$line""_BerryCrust.under100bp" --too-long-output "merged_trimmed_QF/discards/""$line""_BerryCrust.over260bp" -o "merged_trimmed_QF/""$line""_BerryCrust.final.annotated.fastq" "$line""_BerryCrust.collapsed.annotated.fastq"; done < ../../Sample_list.txt

mkdir ../convert_to_fasta
while read line; do  
obiconvert --sanger --fasta-output \
"merged_trimmed_QF/""$line""_BerryCrust.final.annotated.fastq" > "../convert_to_fasta/""$line""_BerryCrust""".fasta""; done < ../../Sample_list.txt

cd ../convert_to_fasta
cat *_BerryCrust.fasta  > BerryCrust.allsamples.fasta
obistat -c sample -a seq_length BerryCrust.allsamples.fasta > sample_stats_BerryCrust.txt 
obisplit -t sample BerryCrust.allsamples.fasta 
while read line; do  
rm "$line"".fasta"; done < ../../Sample_list.txt

#only keep unique sequences in each sample
for i in *_BerryCrust.fasta ; do obiuniq $i >  ${i/.fasta/.unique.fasta}; done

cd ..
mkdir owitools
cd owitools
cp -r ~/Desktop/Metabarcoding_Physalia/ubuntu/owi_tools/ .
chmod 777 *
which RScript 
result on my computer: /usr/local/bin/RScript
copy and paste that into the first line of the owi_script that I want to use
also had to go into RStudio and install.packages("optparse")

#remove chimaeras using vsearch
for i in *_BerryCrust.unique.fasta ; do ./../owitools/owi_obisample2vsearch -i $i; done
for i in *_BerryCrust.unique.vsearch.fasta ; do vsearch --uchime_denovo $i --sizeout --minh 0.90 --nonchimeras ${i/.unique.vsearch.fasta/.nonchimeras.fasta} --chimeras ${i/.unique.vsearch.fasta/.chimeras.fasta} --uchimeout ${i/.unique.vsearch.fasta/.uchimeout.txt}; done


#Concatenating non chimaeras in a single file
cat  *_BerryCrust.nonchimeras.fasta > ../swarm/BerryCrust.nonchimeras.fasta

#Returning to obitools format
../owitools/owi_vsearch2obifasta -i BerryCrust.nonchimeras.fasta

#Dereplicating sequences in MiFish.nonchimeras.fasta file
obiuniq -m sample BerryCrust.nonchimeras.vsearch.fasta > BerryCrust.unique.fasta

#change names of motus to a short more informative identifier
obiannotate --seq-rank BerryCrust.unique.fasta | obiannotate --set-identifier '"BerryCrust_%07d" % seq_rank' > BerryCrust.new.fasta


#obtain a table of abundances
obitab -o BerryCrust.new.fasta >  BerryCrust.new.tab

#back to vsearch format
#../owitools/owi_obifasta2vsearch -i BerryCrust.new.fasta -o BerryCrust.vsearch.fasta #not working right
# I copied a new version of owi_obifasta2vsearch from the original folder and modified it to know that the name was the first 18 characters instead of 14 characters for BerryCrust_%07d
#saved as ../owitools/owi_obifasta2vsearch18

../owitools/owi_obifasta2vsearch18 -i BerryCrust.new.fasta -o BerryCrust.vsearch.fasta



#cluster with swarm :BerryCrust
conda activate eDNA
swarm -d 1 -f -z -t 40 -o BerryCrust.SWARM1nc_output -s BerryCrust.SWARM1nc_stats -w BerryCrust.SWARM1nc_seeds.fasta BerryCrust.vsearch.fasta


swarm -d 2 -z -t 40 -o BerryCrust.SWARM2nc_output -s BerryCrust.SWARM2nc_stats -w BerryCrust.SWARM2nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 3 -z -t 40 -o BerryCrust.SWARM3nc_output -s BerryCrust.SWARM3nc_stats -w BerryCrust.SWARM3nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 4 -z -t 40 -o BerryCrust.SWARM4nc_output -s BerryCrust.SWARM4nc_stats -w BerryCrust.SWARM4nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 5 -z -t 40 -o BerryCrust.SWARM5nc_output -s BerryCrust.SWARM5nc_stats -w BerryCrust.SWARM5nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 6 -z -t 40 -o BerryCrust.SWARM6nc_output -s BerryCrust.SWARM6nc_stats -w BerryCrust.SWARM6nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 7 -z -t 40 -o BerryCrust.SWARM7nc_output -s BerryCrust.SWARM7nc_stats -w BerryCrust.SWARM7nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 8 -z -t 40 -o BerryCrust.SWARM8nc_output -s BerryCrust.SWARM8nc_stats -w BerryCrust.SWARM8nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 9 -z -t 40 -o BerryCrust.SWARM9nc_output -s BerryCrust.SWARM9nc_stats -w BerryCrust.SWARM9nc_seeds.fasta BerryCrust.vsearch.fasta
swarm -d 10 -z -t 40 -o BerryCrust.SWARM10nc_output -s BerryCrust.SWARM10nc_stats -w BerryCrust.SWARM10nc_seeds.fasta BerryCrust.vsearch.fasta
#swarm -d 11 -z -t 40 -o BerryCrust.SWARM11nc_output -s BerryCrust.SWARM11nc_stats -w BerryCrust.SWARM11nc_seeds.fasta BerryCrust.vsearch.fasta
#swarm -d 12 -z -t 40 -o BerryCrust.SWARM12nc_output -s BerryCrust.SWARM12nc_stats -w BerryCrust.SWARM12nc_seeds.fasta BerryCrust.vsearch.fasta

#d values for swarm test y=number of swarms BerryCrust
a <- c(1,2,3,4,5,6,7,8,9,10)
b<- c(8474,7199,4297,
plot(a,b)
atrim <-c(3,4,5,6,7,8,9,10)
btrim<-c(4297,)
plot(atrim,btrim)

../owitools/owi_recount_swarm BerryCrust.SWARM1nc_output BerryCrust.new.tab

#Still!!!!!
Kept only 0 sequences for calculations.


##### Dec 16, 2023: Found comparison of methods for finding intraspecific variation ##### going to try their winning combination:
#method 2B from: https://gitlab.mbb.univ-montp2.fr/edna/exploitation/edna_intra_pipeline_comparison#step22
working in /Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/MF_BC_comparison_pipeline_2b_dec16-22

pre-processing steps:
started with the results of cutadapt and trimmomatic:
    MiFish  /Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/clean_w_cutadapt_dec15-22/trimmed/"$line""_MiFish_R1.trimmed.fastq" 
    BerryCrust  /Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/BC_clean_w_cutadapt_dec15-22/trimmed/"$line""_BerryCrust_R1.trimmed.fastq"

they had been processed yesterday with the following commands:

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

#ok, let's begin with this new method!#

conda activate obitools

#BerryCrust:

while read line; do
illuminapairedend --score-min=40 -r "/Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/BC_clean_w_cutadapt_dec15-22/trimmed/""$line""_BerryCrust_R1.trimmed.fastq" "/Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/BC_clean_w_cutadapt_dec15-22/trimmed/""$line""_BerryCrust_R2.trimmed.fastq" > "$line""_BerryCrust".fastq; done < ../Sample_list.txt
# a new .fastq file is created, it contains the sequences after the merging of forward and reverse strands
# alignments which have a quality score higher than 40 (-- score-min=40) are merged and annotated "aligned", while alignemnts with a lower quality score are concatenated and annotated "joined"

#MiFish
while read line; do
illuminapairedend --score-min=40 -r "/Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/clean_w_cutadapt_dec15-22/trimmed/""$line""_MiFish_R1.trimmed.fastq" "/Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/clean_w_cutadapt_dec15-22/trimmed/""$line""_MiFish_R2.trimmed.fastq" > "$line""_MiFish".fastq; done < ../Sample_list.txt
**********
#BerryCrust
while read line; do
obigrep -p 'mode!="joined"' "$line""_BerryCrust.fastq" > "$line""_BerryCrust.ali.fastq"; done < ../Sample_list.txt
# -p requires a python expression
# python creates a new dataset (.ali.fastq) which only contains the sequences annotated "aligned"

#MiFish
while read line; do
obigrep -p 'mode!="joined"' "$line""_MiFish.fastq" > "$line""_MiFish.ali.fastq"; done < ../Sample_list.txt
# -p requires a python expression
# python creates a new dataset (.ali.fastq) which only contains the sequences annotated "aligned"


#since they're already demultiplexed, I'm going to substitute obiannotate for ngsfilter

mkdir ready_for_dada2

#BerryCrust
for file in *_BerryCrust.ali.fastq; do
sample=${file#*dir/}
sample=${sample%__BerryCrust.ali.fastq}
obiannotate -S sample:${sample} -S assay:BerryCrust --length ${sample}_BerryCrust.ali.fastq > ready_for_dada2/${sample}_BerryCrust.ali.annotated.fastq
done

#MiFish
for file in *_MiFish.ali.fastq; do
sample=${file#*dir/}
sample=${sample%__MiFish.ali.fastq}
obiannotate -S sample:${sample} -S assay:MiFish --length ${sample}_MiFish.ali.fastq > ready_for_dada2/${sample}_MiFish.ali.annotated.fastq
done

made yet another Rscript and stored it in /Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/MF_BC_comparison_pipeline_2b_dec16-22

saved as Mace-Manel-NC_Sharks_comparison_2b.R

output of this file is: PipelineB_", sample.names[k], ".fasta

obigrep -p 'count>=10' Aquarium_2.fasta > Aquarium_2.grep.fasta
# "-p 'count>=10'" option eliminates sequences with an abundance inferior to 10

#they then read this fasta back into dada2 for chimera-removal

#then they don't talk about which method they use for taxonomic assignment.  
rdp? dad/dads? ecotag?
I guess I can try all 3 with the output of NC_Sharks_dada2_dec-16-22.R
stored in BC_dada2_dec-16-22

#while waiting for illumina-paired-end to finish

#rdp has already been done and is waiting for metabar analysis and MiFish to catch up.
#dad/dads dbs have already been prepared and can easily be added
#ecotag is a bit harder.  Need to download the EMBL database (already underway on panthera) and process it with ecopcr
    #look into what format the swarm.seeds.fasta file that worked with ecotag for the fish was.compared to the output of dada2.
    
#also found this OTHER comparison of metabarcoding pipeline steps that includes taxonomic classifiers!#
https://www.vigilife.org/wp-content/uploads/mathon-et-al-2021-molecolres.pdf
#from https://github.com/lmathon/eDNA--benchmark_pipelines/blob/master/benchmark_real_dataset/06_assignation/Assign_vsearch.sh

#they recommend VSEARCH –usearch_global for the taxonomic assignment.

# Taxonomic assignation
all_sample_sequences_vsearch_tag="${all_sample_sequences_RC/.fasta/.tag.fasta}"
$vsearch --usearch_global $all_sample_sequences_RC --db $refdb_dir --qmask none --dbmask none --notrunclabels --id 0.98 --top_hits_only --threads 16 --fasta_width 0 --maxaccepts 20 --maxrejects 20 --minseqlength 1 --maxhits 20 --query_cov 0.6 --blast6out $all_sample_sequences_vsearch_tag --dbmatched $main_dir/db_matched.fasta --matched $main_dir/query_matched.fasta
## Create final table
### preformat
all_sample_sequences_vsearch_preformat="${all_sample_sequences_vsearch_tag/.fasta/.preformat.fasta}"
tr "\t" ";" < $all_sample_sequences_vsearch_tag > $all_sample_sequences_vsearch_preformat
python3 /benchmark_real_dataset/optimized_pipeline/vsearch2obitab.py -a $all_sample_sequences_vsearch_preformat -o $fin_dir/opt_pipeline.csv

#also from https://github.com/lmathon/eDNA--benchmark_pipelines/blob/master/benchmark_real_dataset/06_assignation/Assign_vsearch.sh

test dataset
vsearch --usearch_global $all_sample_sequences_uniqid --db $refdb_dir --qmask none --dbmask none --notrunclabels --id 0.98 --top_hits_only --threads 16 --fasta_width 0 --maxaccepts 20 --maxrejects 20 --minseqlength 20 --maxhits 20 --query_cov 0.6 --blast6out $all_sample_sequences_vsearch_tag --dbmatched $main_dir/db_matched.fasta --matched $main_dir/query_matched.fasta

python3 benchmark_real_dataset/06_assignation/vsearch2obitab.py -a $all_sample_sequences_vsearch_tag -o $fin_dir/"$step".csv
#from vsearch documentation
vsearch --usearch_global queries.fas --db references.fas --id 0.8 --iddef 1 --alnout results.aln


###### Before I really got into trying the comparison pipeline or the vsearch taxonomic assignment of the dada2 reads... BerryCrust worked with the obitools3 method!!!!! ######
halleluiah!

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

the rest is in FINAL_Obitools3 processing of NC_Sharks project.txt (beginning to end) 
stored on the NAS drive in the panthera_data folder
and also in ~/Desktop/Savannah/NC_Sharks_data/Cutadapt
and  ~Desktop/Savannah/ (latest version with excel and metabar notes)

12-18

after FINAL_Obitools3 processing, I imported the fish_results.fasta and crustaceans_results.fasta to my computer in: /Users/Eldridge/Desktop/Savannah/following_obitools3_method

#Change identifiers by a short index
obiannotate --seq-rank fish_results.fasta | obiannotate --set-identifier '"MiFish_%07d" % seq_rank' > MiFish_named.fasta
obiannotate --seq-rank crustaceans_results.fasta | obiannotate --set-identifier '"BerryCrust_%07d" % seq_rank' > BerryCrust_named.fasta

#obtain a table of abundances
obitab -o MiFish_named.fasta >  MiFish_named.tab
obitab -o BerryCrust_named.fasta >  BerryCrust_named.tab

#12-21-22 Try lulu directly after all obitools3 pipeline.

Used BerryCrust_named.tab to make the otu table for lulu by removing all of the extra columns that lulu doesn't use (leaving OTU_ID (labeled it OBI3_ID), and the abundance table)and changing all NAs in the table to 0, and removing MERGED_sample: from the front of every sample name in the table.  Did this in excel.
saved it as BerryCrust_named_tab_lulu_formatted.txt
did the lulu curation in post-dada2-obi3_lulu.R file.
saved resulting file as: lulu_obi3_BerryCrust_named_results_otutable.tsv

#now I have to get the taxonomic information back into the otutable (owi_combine?), or do it in MetabaR or phyloseq

#added the sample. to the beginning ofthe lulu output file samplenames in excel and saved as:
lulu_obi3_BerryCrust_named_results_otutable_forowicombine.txt
#it also wants and ecotag taxonomy labeled csv possibly I can make one from BerryCrust_named.tab in excel.
deleted the otutable from BerryCrust_named.tab and saved as BerryCrust_named_tab_taxonomyforowicombine.csv

ok, here goes!

./owitools/owi_combine -i BerryCrust_named_tab_taxonomyforowicombine.csv -a lulu_obi3_BerryCrust_named_results_otutable_forowicombine.txt -o BerryCrust_named_lulu_with_tax.txt

Ecotag database read including 2875 total MOTUs.
Reading abundance database...
Abundances database read including 839 total MOTUs and  samples.
Error in fix.by(by.x, x) : 'by' must specify a uniquely valid column
Calls: merge -> merge.data.frame -> fix.by
Execution halted

renamed OBI3_ID to id since the code says fix.by "id"
#didn't solve the error. checking if the csv and txt comma and tab delimiters are causing the problem:
saved the BerryCrust_named_tab_taxonomyforowicombine file as txt tab delimited
./owitools/owi_combine -i BerryCrust_named_tab_taxonomyforowicombine.txt -a lulu_obi3_BerryCrust_named_results_otutable_forowicombine.txt -o BerryCrust_named_lulu_with_tax.txt --sep="\t"
Reading ecotag database...
Ecotag database read including 2875 total MOTUs.
Reading abundance database...
Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  : 
  line 1 did not have 101 elements
Calls: read.table -> scan
Execution halted

oh well, I'll try to do it myself in R instead. 
OR, I could transform it into a fasta and put it back through ecotag!  

let's clean up the owicombine input files.
rm *forowicombine*

#writing lulu_obi3_BerryCrust_named_results_otutable.tsv to fasta format for ecotag
OR I could put it in phyloseq as-is...?!

trying that with Phyloseq_obi3_lulu.R
I found that the BerryCrust_named.fasta and .tab files only have ids to genus!  how did that happen?  I will re-export from obi3


12-22-21 made a new directory in following_obitools3_method called results_obi3-12-21-22 and copied
the crustaceans_results.fasta and .tsv files from the re-export to panthera_data yesterday.  
checked that they contain annotations with 
grep -c 'Callinectes similis' crustaceans_results.tsv
62

tried to import them into R for lulu curation but it says duplicate rownames aren't allowed.  I'm sure it's just a problem with how R is reading the names
I'm going to try to rename them but keep the annotations.








##########################################################################
#12-18-22
#exported the good_sequences for the fish and crustacean DMS views 
#copied them to my computer to use owitools and vsearch to remove chimeras then re-imported as good_sequences_nochim

###############   putting human-readable unique identifiers and removing chimaeras with vsearch
obi export --fasta-output fish/good_sequences -o MiFish_aligned.fasta
obi export --fasta-output crustaceans/good_sequences -o BerryCrust_aligned.fasta

copy to my computer
#split by sample
conda activate obitools

obisplit -t sample -p MiFish_ MiFish_aligned.fasta
obisplit -t sample -p BerryCrust_ BerryCrust_aligned.fasta

for i in BerryCrust_*.fasta ; do obiuniq $i >  ${i/.fasta/.unique.fasta}; done
for i in MiFish_*.fasta ; do obiuniq $i >  ${i/.fasta/.unique.fasta}; done


change to vsearch format to remove chimeras on my computer using owi_tools

for i in BerryCrust_*.unique.fasta ; do ./owitools/owi_obisample2vsearch -i $i; done
for i in MiFish_*.unique.fasta ; do ./owitools/owi_obisample2vsearch -i $i; done

conda activate eDNA
#version = vsearch v2.15.2_macos_x86_64
for i in BerryCrust_*.unique.vsearch.fasta ; do vsearch --uchime_denovo $i --sizeout --minh 0.90 --nonchimeras ${i/.unique.vsearch.fasta/.nonchimeras.fasta} --chimeras ${i/unique.vsearch.fasta/.chimeras.fasta} --uchimeout ${i/unique.vsearch.fasta/.uchimeout.txt}; done
for i in MiFish_*.unique.vsearch.fasta ; do vsearch --uchime_denovo $i --sizeout --minh 0.90 --nonchimeras ${i/.unique.vsearch.fasta/.nonchimeras.fasta} --chimeras ${i/unique.vsearch.fasta/.chimeras.fasta} --uchimeout ${i/unique.vsearch.fasta/.uchimeout.txt}; done



#Concatenating non chimaeras in a single file
cat  BerryCrust_*.nonchimeras.fasta > outputs/BerryCrust.nonchimeras.fasta
cat  MiFish_*.nonchimeras.fasta > outputs/MiFish.nonchimeras.fasta

cd outputs/
#Returning to obitools format
../owitools/owi_vsearch2obifasta -i BerryCrust.nonchimeras.fasta
../owitools/owi_vsearch2obifasta -i MiFish.nonchimeras.fasta

obiannotate -S assay:MiFish  MiFish.nonchimeras.vsearch.fasta > MiFish.nonchimeras.annotated.fasta
obiannotate -S assay:BerryCrust  BerryCrust.nonchimeras.vsearch.fasta > BerryCrust.nonchimeras.annotated.fasta


#Dereplicating sequences/motus-to-be in MiFish.nonchimeras.fasta file
obiuniq -m sample BerryCrust.nonchimeras.annotated.fasta > BerryCrust.unique.fasta
obiuniq -m sample MiFish.nonchimeras.annotated.fasta > MiFish.unique.fasta



#change names of motus to a short more informative identifier
obiannotate --seq-rank BerryCrust.unique.fasta | obiannotate --set-identifier '"BerryCrust_%07d" % seq_rank' > BerryCrust.new.fasta
obiannotate --seq-rank MiFish.unique.fasta | obiannotate --set-identifier '"MiFish_%07d" % seq_rank' > MiFish.new.fasta


#obtain a table of abundances
obitab -o BerryCrust.new.fasta >  BerryCrust.new.tab
obitab -o MiFish.new.fasta >  MiFish.new.tab

cat *.new.fasta > NC-Sharks.new.fasta
obitab -o NC-Sharks.new.fasta >NC-Sharks.new.tab
##################################################



formatted the desired 4 tables in Excel.  #I used MiFish.swarm4.LULU.curated.csv, and NCSharks_MiFish_results.fasta (from obi3 tax_assigned)# and the sample sheet for sequencing and the spreadsheets of the sample info(both sampling and lab)

I created the table PCRs by looking at the soil_euk example data and filling in the information from the sample sheet for sequencing
I created a new column in pcrs for sample_id and made it Shark# with the number corresponding to the NC fecal swab number when available and made new Shark#s for unmatched stomachs.  I put these in the Stomach Content Data table that Savannah sent me and put it in the shared google drive
I created another new column in pcrs for material and named them either fecal, stomach, or negative
I added columns for type sample or control and control_type the controls in this experiment were all PCR negatives
I changed IJKMNP of plate 2 (which were the names of the adapters, to ABCDEF)

I sorted Motus and Reads by Motu# order then copied and pasted reads (number of reads of each motu in each sample from LULU curated swarm4 file)
into a new excel sheet but used paste-special, transpose to get it so that the samples were columns and motus were rows.
no first cell label on either of these is needed (don't use one)

Importing caused some issues if the files weren't copied and pasted into new sheets because of empty row or column names/cells at the ends of the data.

Eventually, it worked!

Now, continuing with the tutorial! https://metabarfactory.github.io/metabaR/articles/metabaRF-vignette.html#tutorial-with-the-soil_euk-dataset

12/18/22 got to the end of the metabar analysis pipeline for both but I don't like how the sequences are not clustered by taxa there are so many callinectes columns for example!

12-19-22 I tried to make the dada2 ASV method to obi3 taxonomic classification work again.  I used the R script: /Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/BC_dada2_dec-16-22/NC_Sharks_dada2_dec-16-22.R

then I added a bit of code to extract the unique ASVs from the mergers (derep, dada, merge) made using the latest cleaning (~/Desktop/Savannah/NC_Sharks_data/Cutadapt/BC_clean_w_cutadapt_dec15-22/trimmed/)
and put them in a file for each sample, because that's the only way I could see to get the ASV sequences out with sample information intact.  The only part it's missing is extract bimera denovo

These are now stored in: /Users/Eldridge/Desktop/Savannah/dada2_to_obi3/individual_dada2_asv_fastas/BerryCrust_dec-19-22/

next, I need to , possibly remove chimeras with vsearch, annotate them with sample and assay again, then send to obi3 for taxonomic assignment

conda activate eDNA
for i in *_BerryCrust_ASVs.fasta ; do vsearch --uchime_denovo $i --sizeout --minh 0.90 --nonchimeras ${i/.fasta/.nonchimeras.fasta} --chimeras ${i/.fasta/.chimeras.fasta} --uchimeout ${i/.fasta/.uchimeout.txt}; done

conda activate obitools
for file in *.nonchimeras.fasta; do
sample=${file#*dir/}
sample=${sample%_BerryCrust_ASVs.nonchimeras.fasta}
obiannotate -S sample:${sample} -S assay:BerryCrust --length ${sample}_BerryCrust_ASVs.nonchimeras.fasta > ${sample}_BerryCrust_ASVs.nonchimeras.annotated.fasta
done

cat *_BerryCrust_ASVs.nonchimeras.annotated.fasta >BerryCrust_ASVs.nonchimeras.annotated.fasta
obiannotate --seq-rank BerryCrust_ASVs.nonchimeras.annotated.fasta | obiannotate --set-identifier '"BerryCrust_%08d" % seq_rank' > BerryCrust_ASVs.nonchimeras.annotated.new.fasta

#another crossroads...

#not sure this obitab will ever be needed.
obitab -o BerryCrust_ASVs.nonchimeras.annotated.new.fasta >  BerryCrust_ASVs.nonchimeras.annotated.new.tab


####import into obi3 do obi unique, ecotag, stats
imports: /Users/Eldridge/Desktop/Savannah/dada2_to_obi3/individual_dada2_asv_fastas/BerryCrust_dec-19-22
notes on panthera
exports: /Users/Eldridge/Desktop/Savannah/dada2_to_obi3/results
####


create an asv table for LULU
#vsearch --search_exact reads.fasta #everything before uniques to fasta from dada2 -db asvs.fasta #output of obi3 asvs -otutabout output.tsv
vsearch --search_exact reads.fasta -db asvs.fasta -otutabout output.tsv


exported the merged (pre-uniques) sequences by sample into /Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/BC_dada2_dec-16-22
using the R script: /Users/Eldridge/Desktop/Savannah/dada2_to_obi3/results/NC_Sharks_dada2_to_obi3.R
copied and pasted the merged (pre-uniques) fastas from /Users/Eldridge/Desktop/Savannah/NC_Sharks_data/Cutadapt/BC_dada2_dec-16-22 to /Users/Eldridge/Desktop/Savannah/dada2_to_obi3/individual_dada2_asv_fastas/BerryCrust_dada_pre-uniques

#annotate them all with sample=
for file in *.dada2out.merged.fasta; do
sample=${file#*dir/}
sample=${sample%.dada2out.merged.fasta}
obiannotate -S sample:${sample} ${sample}.dada2out.merged.fasta > ${sample}.dada2out.merged.annotated.fasta
done

switch to vsearch format:
for i in *.dada2out.merged.annotated.fasta ; do ../../owitools/owi_obisample2vsearch -i $i; done

cat *dada2out.merged.annotated.vsearch.fasta > merged_dada2_reads.fasta 
cd /Users/Eldridge/Desktop/Savannah/dada2_to_obi3/results
vsearch --search_exact /Users/Eldridge/Desktop/Savannah/dada2_to_obi3/individual_dada2_asv_fastas/BerryCrust_dada_pre-uniques/merged_dada2_reads.fasta  -db crustaceans_ASV_results.fasta -otutabout crustaceans_ASV_otutable.tsv
#re-do with vsearch --usearch_global merged.fasta --db otus.fasta --id 0.9 --otutabout otutab.txt because the lulu curation didn't do anything.
vsearch --usearch_global /Users/Eldridge/Desktop/Savannah/dada2_to_obi3/individual_dada2_asv_fastas/BerryCrust_dada_pre-uniques/merged_dada2_reads.fasta --db crustaceans_ASV_results.fasta --id 0.9 --otutabout crustaceans_ASV_otutable.tsv

#create matchlist for LULU
vsearch --usearch_global crustaceans_ASV_results.fasta --db crustaceans_ASV_results.fasta --self --id .84 --iddef 1 --userout crustaceans_ASV_results.match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10

ran lulu on these files.  post-dada2-obi3_lulu.R. 

#next I will build the metabarlist using the results of lulu as the table

use lulu_curated_crustaceans_ASV_otutable.tsv and crustaceans_ASV_results.csv to make the files for MetabaR
##!!!!!
#this doesn't seem right... I probably shouldn't have done uniques twice... maybe I should try this without doing uniques in the obi3pipeline.


#did this same thing for the obi3 pure(non-dada2)files

cd /Users/Eldridge/Desktop/Savannah/following_obitools3_method
check both files for annotations and formatting 

sed -i -e 's/ //g' crustaceans_merged_reads.fasta

vsearch --usearch_global crustaceans_merged_reads.fasta --db BerryCrust_named.fasta --id 0.9 --otutabout obi3_crustaceans_results_otutable.tsv

#create a matchlist for LULU (not the one from the obi 3 align command)
vsearch --usearch_global BerryCrust_named.fasta --db BerryCrust_named.fasta --self --id .84 --iddef 1 --userout crustaceans_results_BerryCrust_named_matchlist.tsv -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10

ran lulu on these files.  post-dada2-obi3_lulu.R (second section)

result:lulu_obi3_crustaceans_results_otutable.tsv use this and BerryCrust_named.fasta to make the MetabaR tables

St9 was lost before the read-in of obi3_crustaceans_results_otutable.tsv and before BerryCrust_named.fasta 

and I don't see a file called St9 in BC_clean_w_cutadapt_dec15-22 folder even
deleted it from PCRs for lulu reads metabar file

same issue with the ASVs not wanting to go into metabar
I can just try visualizations fishfarm tomorrow or phyloseq
#couldn't get themtomatch up,wenttobedat 2AM

#could also try taxonomic assignment with vsearch on the custom database or the EMBL ecopcr output 

#I can try the owitools combine after lulu script but it needs a taxonomy file.

#alternatively to these packages, I can import the obi3 results and the lulu counts into R as dataframes
and just curate them with a custom R script, maybe go back to the class work from Physalia Metabarcoding workshop
and do the visualizations from the fishfarm R script.  They curated in excel, deleting nontarget ASVs etc.
except without the dada2 input (which I now think I messed up by doing obi uniq in obi3 because it got rid of the original count information which lulu probably needed)

I guess I can call the obi clean output ASVs...?

####### ok, the best results (I think,so far have been the obitools3 results- dada2_to_obi3 is a close second.)

I found that the BerryCrust_named.fasta and .tab files only have ids to genus!  how did that happen?  I will re-export from obi3
#######################

12-22-21 made a new directory in following_obitools3_method called results_obi3-12-21-22 and copied
the crustaceans_results.fasta and .tsv files from the re-export to panthera_data yesterday.  
checked that they contain annotations with 
grep -c 'Callinectes similis' crustaceans_results.tsv
62

tried to import them into R for lulu curation but it says duplicate rownames aren't allowed.  I'm sure it's just a problem with how R is reading the names
I'm going to try to rename them but keep the annotations.

I'm going to try lulu curation directly on the outputs of obitools3 end results and see if that can get down to one ortwo sequences assigned to Callinectes instead of hundreds

obiannotate --seq-rank crustaceans_results.fasta | obiannotate --set-identifier '"BerryCrust_%07d" % seq_rank' > BerryCrust_named.fasta
grep -c 'Callinectes similis' BerryCrust_named.fasta 
62
#otutable for lulu
obitab -o -d -n 0 BerryCrust_named.fasta >BerryCrust_named.tab
grep -c 'Callinectes similis' BerryCrust_named.tab 
62
make a new tab file with only the info Lulu wants.
#I put the file into excel and added a row to count the field numbers to include
cut -f1,7-164 BerryCrust_named.tab > BerryCrust_named_tab_LULU.txt && sed -i -e 's/MERGED_sample://g' BerryCrust_named_tab_LULU.txt


#matchlist for lulu
#First produce a blastdatabase with the OTUs
obiannotate -C BerryCrust_named.fasta >BerryCrust_named_cleared_tags.fasta
conda activate blast
makeblastdb -in BerryCrust_named_cleared_tags.fasta -parse_seqids -dbtype nucl
# Then blast the OTUs against the database
blastn -db BerryCrust_named_cleared_tags.fasta -outfmt '6 qseqid sseqid pident' -out BerryCrust_named_matchlist.txt -qcov_hsp_perc 80 -perc_identity 84 -query BerryCrust_named_cleared_tags.fasta


used BerryCrust_named_matchlist.txt and BerryCrust_named_tab_LULU.txt in lulu
Phyloseq_obi3_lulu.R file.
lulu with default parameters kept 
654 motus

minimum_ratio_type	
sets whether a potential error must have lower abundance than the parent in all samples min (default), or if an error just needs to have lower abundance on average avg. Choosing lower abundance on average over globally lower abundance will greatly increase the number of designated errors. This option was introduced to make it possible to account for non-sufficiently clustered intraspecific variation, but is not generally recommended, as it will also increase the potential of cluster well-separated, but co-occuring, sequence similar species.
minimum_ratio_type = "avg" result =406 retained motus 

I'll try the default parameters first.
I need to get rid of the pcr contaminants and tag jumps before gettting too far into the phyloseq visualizations
metabar!
worked eventually as soon as I realized lulu was substituting . instead of - to the St#-St# pcr sample names.

#######
Now do the same thing (pure obi3 results, no matchback) for MiFish results.

#checked the input files:
following_obitools3_method Eldridge$ grep -c 'Micropogonias furnieri' fish_results.fasta
13
(obitools) Eldridge-Wiselys-2022-MacbookPro:following_obitools3_method Eldridge$ grep -c 'Micropogonias furnieri' MiFish_named.fasta 
13
(obitools) Eldridge-Wiselys-2022-MacbookPro:following_obitools3_method Eldridge$ grep -c 'Micropogonias furnieri' MiFish_named.tab
13
(obitools) Eldridge-Wiselys-2022-MacbookPro:following_obitools3_method Eldridge$ grep -c 'Micropogonias furnieri' MiFish_named.fasta
13
make a new tab file with only the info Lulu wants.
#I put the file into excel and added a row to count the field numbers to include MiFish_named_tab.xlsx
cut -f1,8-165 MiFish_named.tab > MiFish_named_tab_LULU.txt && sed -i -e 's/MERGED_sample://g' MiFish_named_tab_LULU.txt
sed -i -e 's/NA/0/g' MiFish_named_tab_LULU.txt

#so this doesn't get even more confusing...
mkdir MiFish_pureobi3_lulu_metabar
cp MiFish_named* MiFish_pureobi3_lulu_metabar/
cd MiFish_pureobi3_lulu_metabar/

#matchlist for lulu
#First produce a blastdatabase with the OTUs
obiannotate -C MiFish_named.fasta >MiFish_named_cleared_tags.fasta
conda activate blast
makeblastdb -in MiFish_named_cleared_tags.fasta -parse_seqids -dbtype nucl
# Then blast the OTUs against the database
blastn -db MiFish_named_cleared_tags.fasta -outfmt '6 qseqid sseqid pident' -out MiFish_named_matchlist.txt -qcov_hsp_perc 80 -perc_identity 84 -query MiFish_named_cleared_tags.fasta


used MiFish_named_matchlist.txt and MiFish_named_tab_LULU.txt in lulu
made a new MiFish_Phyloseq_obi3_lulu.R file. Based on the Phyloseq_obi3_lulu.R file that I used to process the MiFish reads successfully.

lulu with default parameters kept 

12/29/22
DONE!










