#!/bin/bash
echo "Script start: download and initial sequencing read quality control" 
date
#Navigate to analysis folder in shared project directory
cd medbioinfo_folder/lina/MedBioinfo/analyses
#Save the run numbers of chosen data
sqlite3 -batch -noheader -csv ../../../pascal/central_database/sample_collab.db "SELECT run_accession FROM sample_annot LEFT JOIN sample2bioinformatician ON sample_annot.patient_code=sample2bioinformatician.patient_code WHERE username = 'linastep';" > linastep_run_accessions.txt #results in text file with one column (not 2 as said in the assignment)
#Create fastq subdirectory in data folder
mkdir ../data/sra_fastq
#Load module to download the files
module load sra-tools
#Download the files with accessions
cat linastep_run_accessions.txt | srun --cpus-per-task=1 --time=01:00:00 xargs fastq-dump --split-3 --gzip --readids --outdir ../data/sra_fastq/ --disable-multithreading
#count how many reads were found per FASTQ file
find ../data/sra_fastq/ -name "*.fastq.gz" | xargs zgrep -c "^@"
#compare number of acquired reads with seqkit
module load seqkit/
srun --cpus-per-task=4 --time=01:00:00 seqkit stats -j 4 ../data/sra_fastq/*.fastq.gz > linastep_seqkit_stats.txt
sqlite3 -batch ../../../pascal/central_database/sample_collab.db "SELECT * FROM sample_annot LEFT JOIN sample2bioinformatician ON sample_annot.patient_code=sample2bioinformatician.patient_code WHERE username = 'linastep';" #after manual comparison number of reads correspond; the reads are trimmed as normally Illumina reads are 50-300 bp length and the reads that we have are 35-151 bp length. 
#check if FASTQ files have been de-replicated. 

echo "Finished"
date
