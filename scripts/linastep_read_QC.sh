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


echo "Finished"
date
