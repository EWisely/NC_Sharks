# NC_Sharks
Code used in NC_Sharks project [^1] for shark diet metabarcoding from fecal swabs and stomach contents. The laboratory protocol for creating these data can be found on protocols.io at DOI: dx.doi.org/10.17504/protocols.io.ewov1qxokgr2/v1.[^2] This github repo covers the analysis steps and contains 6 scripts to be used in order:
1) 01_demultiplex_and_clean.sh demultiplexes the two metabarcoding primers we used: MiFish[^3] and Berry's crustacean primer (BerryCrust in our code)[^4] and removes the adapters from the raw data from the MiSeq machine using Adapterama iNext primers [^5] and Illumina overhang adapters[^6].
2) 02_Obitools3_Metabarcoding.sh takes the cleaned, annotated R1 and R2 fastq files for both the crustaceans and fishes and metabarcodes them using Obitools3[^7] for merging, initial filtering, and cleaning with obiclean to determine the "head" sequences (MOTUs) to be taxonomically identified, then compares these obicleaned "head" sequences to the EMBL databases for invertebrates and vertebrates respectively using the ecopcr[^8] and ecotag functions.  The outputs are crustaceans_results.fasta and fish_results.fasta
3) 03_prep_for_LULU.sh uses the output fasta files of taxonomically identified MOTUs from Obitools3 as input, and uses Obitools (version 2), and BLAST[^9],[^10] to create the input files for the R[^11] scripts.  The output files from this script are: MiFish_named_matchlist.txt, MiFish_named_tab_LULU.txt, BerryCrust_named_matchlist.txt, and BerryCrust_named_tab_LULU.txt
4) 04_lulu_metabaR_BerryCrust_NC_Sharks.R takes as input BerryCrust_named_matchlist.txt, BerryCrust_named_tab_LULU.txt, containing the metabarcoding results, and MetabaR_formatted_samples_NCSharks.txt, and MetabaR_formatted_pcrs_NCSharks_BerryCrust.txt containing metadata collected during shark captures, and during the library building process respectively.  Using LULU[^12],[^13], taxonomizr [^14], MetabaR[^15],  and Phyloseq[^16] R packages, it curates and cleans the metabarcoding MOTUs to remove spurious signal using data from negative controls, replicates, sequencing depth, and sequence similarity scores.  The results of this script are the final crustacean reads from the fecal and stomach samples.
5) 




[^1]: Ryburn et al. (in preparation).
[^2]: Wisely, Eldridge. 2023. “Metabarcoding Fecal Swabs or Stomach Contents for Fish and Crustaceans Using 2-PCR Protocol and Illumina MiSeq.” Protocols.Io. 2023. https://www.protocols.io/private/AF851F699D1B11EE9FCB0A58A9FEAC02.
[^3]: Miya, M., Y. Sato, T. Fukunaga, T. Sado, J. Y. Poulsen, K. Sato, T. Minamoto, et al. 2015. “MiFish, a Set of Universal PCR Primers for Metabarcoding Environmental DNA from Fishes: Detection of More Than 230 Subtropical Marine Species.” Royal Society Open Science 2 (7): 150088. https://doi.org/10.1098/rsos.150088.
[^4]: Berry, Tina E., Benjamin J. Saunders, Megan L. Coghlan, Michael Stat, Simon Jarman, Anthony J. Richardson, Claire H. Davies, Oliver Berry, Euan S. Harvey, and Michael Bunce. 2019. “Marine Environmental DNA Biomonitoring Reveals Seasonal Patterns in Biodiversity and Identifies Ecosystem Responses to Anomalous Climatic Events.” PLoS Genetics 15 (2). https://doi.org/10.1371/journal.pgen.1007943.
[^5]: Glenn, Travis C., Todd W. Pierson, Natalia J. Bayona-Vásquez, Troy J. Kieran, Sandra L. Hoffberg, Jesse C. Thomas, Daniel E. Lefever, et al. 2019. “Adapterama II: Universal Amplicon Sequencing on Illumina Platforms (TaggiMatrix).” PeerJ 2019 (10). https://doi.org/10.7717/peerj.7786.
[^6]:  “16S Metagenomic Sequencing Library Preparation.” 2013. Illumina. https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf.
[^7]:Boyer F, Mercier C, Bonin A, Le Bras Y, Taberlet P, Coissac E: OBITools: a Unix-inspired software package for DNA metabarcoding. Mol Ecol Resour, 2016: 176-182.
[^8]: Ficetola, Gentile Francesco, Eric Coissac, Stéphanie Zundel, Tiayyba Riaz, Wasim Shehzad, Julien Bessière, Pierre Taberlet, and François Pompanon. 2010. “An In Silico Approach for the Evaluation of DNA Barcodes.” BMC Genomics 11 (July): 434. 
[^9]: Altschul, S. F., W. Gish, W. Miller, E. W. Myers, and D. J. Lipman. 1990. “Basic Local Alignment Search Tool.” Journal of Molecular Biology 215 (3): 403–10. https://doi.org/10.1016/S0022-2836(05)80360-2.
[^10]: Camacho, Christiam, George Coulouris, Vahram Avagyan, Ning Ma, Jason Papadopoulos, Kevin Bealer, and Thomas L. Madden. 2009. “BLAST+: Architecture and Applications.” BMC Bioinformatics 10 (1): 421. https://doi.org/10.1186/1471-2105-10-421.
[^11]: R Core Team (2022). R: A language and environment for statistical computing. R
  Foundation for Statistical Computing, Vienna, Austria. URL
  https://www.R-project.org/.
[^12]: Guldberg Frøslev T (2022). _lulu: Post Clustering Curation of Amplicon Data for
  Reliable Biodiversity Metrics_. R package version 0.1.0.
[^13]: Frøslev, Tobias Guldberg, Rasmus Kjøller, Hans Henrik Bruun, Rasmus Ejrnæs, Ane Kirstine Brunbjerg, Carlotta Pietroni, and Anders Johannes Hansen. 2017. “Algorithm for Post-Clustering Curation of DNA Amplicon Data Yields Reliable Biodiversity Estimates.” Nature Communications 8 (1). https://doi.org/10.1038/s41467-017-01312-x.
[^14]: Sherrill-Mix S (2023). _taxonomizr: Functions to Work with NCBI Accessions and
  Taxonomy_. R package version 0.10.2,
  <https://CRAN.R-project.org/package=taxonomizr>.
[^15]: Lucie Zinger, Clément Lionnet, Anne-Sophie Benoiston, Julian Donald, Céline
  Mercier, Frédéric Boyer (2021). metabaR : an R package for the evaluation and
  improvement of DNA metabarcoding data quality. Methods in Ecology and Evolution;
  doi: https://doi.org/10.1111/2041-210X.13552
[^16]: McMurdie, Paul J., and Susan Holmes. 2013. “Phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data.” Edited by Michael Watson. PLoS ONE 8 (4): e61217. https://doi.org/10.1371/journal.pone.0061217.

  
