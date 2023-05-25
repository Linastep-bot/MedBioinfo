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

#check if FASTQ files have been de-replicated. Trying it on only one of the FASTQ files as it is safe to assume that if there are duplicates in one it will be in the whole study.
zcat ../data/sra_fastq/ERR6913147_1.fastq.gz | srun --cpus-per-task=4 --time=00:30:00 seqkit rmdup -j 4 -s -o ../data/sra_fastq/clean_ERR6913147_1.fastq.gz #returned 109446 duplicates removed, therefore we can claim that reads are not de-replicated. Given that we want to see the qualitative increase in some of the reads compared to qualitatitve number of different reads, non de-replicated file types are better to be used. 

#check if adapters have been removed. Trying it on only one of the FASTQ files as it is safe to assume that if there are duplicates in one it will be in the whole study.
zcat ../data/sra_fastq/ERR6913147_1.fastq.gz | srun --cpus-per-task=4 --time=00:30:00 seqkit locate -j 4 -p AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -M
zcat ../data/sra_fastq/ERR6913147_1.fastq.gz | srun --cpus-per-task=4 --time=00:30:00 seqkit locate -j 4 -p AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -M #neither of full adapters returned any matches
zcat ../data/sra_fastq/ERR6913147_1.fastq.gz | srun --cpus-per-task=4 --time=00:30:00 seqkit locate -j 4 -p AGATCGGAAGAGCACAC -M #returned several matching sites but none started at the beggining of the sequence (not an adapter, just natural sequence of the strand)
zcat ../data/sra_fastq/ERR6913147_1.fastq.gz | srun --cpus-per-task=4 --time=00:30:00 seqkit locate -j 4 -p AGATCGGAAGAGCGTCG -M #tried approx half-sites of the read2 adapter, no match due to the fact that adaptor is for - strands. 

#Quality control of FASTQ files
srun --cpus-per-task=2 --time=00:30:00 fastqc -o ./fastqc/ -t 2 --noextract ../data/sra_fastq/ERR6913147_1.fastq.gz ../data/sra_fastq/ERR6913147_2.fastq.gz
#running all files at the same time
srun --cpus-per-task=2 --time=00:30:00 xargs -I{} -a linastep_run_accessions.txt fastqc -o ./fastqc/ -t 2 --noextract ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz
#according to the fastqc files, the low quality sequences and sequences adapters were removed. 

#Merging pair-end reads using flash2
srun --cpus-per-task=2 --time=00:30:00 flash2 --threads=2 --output-directory=../data/merged_pairs --output-prefix=ERR6913147.flash --compress ../data/sra_fastq/ERR6913147_1.fastq.gz ../data/sra_fastq/ERR6913147_2.fastq.gz 2>&1 | tee -a linastep_flash2.log #64% of the reads were merged successfully; program suggest increasing -M parameter to be bigger than 65

#checking the length of the merged reads
srun --cpus-per-task=1 --time=00:30:00 seqkit stats -j 1 ../data/merged_pairs/ERR6913147.flash.extendedFrags.fastq.gz #longest read now is 292 instead of 160


#merging all paired reads 
srun --cpus-per-task=2 --time=00:30:00 xargs -a linastep_run_accessions.txt -n 1 -I{} flash2 --threads=2 --output-directory=../data/merged_pairs --output-prefix={}.flash --compress ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz 2>&1 | tee -a linastep_flash2.log

#Removing PhiX contamination 
mkdir ../data/reference_seqs #making directory to store reference sequences

sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)" #installing the ncbi package in working folder rather than home directory

efetch -db nuccore -id NC_001422 -format fasta > ../data/reference_seqs/PhiX_NC_001422.fna #getting NC_001422 sequence as a test

module load bowtie2 #load bowtie module

#converting reference sequence to bowtie database format 
mkdir ../data/bowtie2_DBs
srun bowtie2-build -f ../data/reference_seqs/PhiX_NC_001422.fna ../data/bowtie2_DBs/PhiX_bowtie2_DB

#checking whether there is PhiX contamination
mkdir bowtie
srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_DBs/PhiX_bowtie2_DB -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz -S bowtie/linastep_merged2PhiX.sam --threads 8 --no-unal 2>&1 | tee bowtie/linastep_bowtie2_PhiX.log #none were found

#aligning the reads to the Covid reference gene
efetch -db nuccore -id NC_045512 -format fasta > ../data/reference_seqs/SC2_NC_045512.fna
srun bowtie2-build -f ../data/reference_seqs/SC2_NC_045512.fna ../data/bowtie2_DBs/SC2_bowtie2_DB
srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_DBs/SC2_bowtie2_DB -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz -S bowtie/linastep_merged2SC2.sam --threads 8 --no-unal 2>&1 | tee bowtie/linastep_bowtie2_SC2.log #3 reads were aligned to the Sars-Cov-2 reference genome

#aligning to most promanent nasal bacterial specie Corynebacterium accolens
efetch -db nuccore -id NR_042139 -format fasta > ../data/reference_seqs/Ca_NR_042139.fna
srun bowtie2-build -f ../data/reference_seqs/Ca_NR_042139.fna ../data/bowtie2_DBs/Ca_bowtie2_DB
srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_DBs/Ca_bowtie2_DB -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz -S bowtie/linastep_merged2SCa.sam --threads 8 --no-unal 2>&1 | tee bowtie/linastep_bowtie2_Ca.log #20288 reads (0.92%) were aligned to microbial reference sequence

#perform a multiQC on all files 
srun multiqc --force --title "linastep sample sub-set" ../data/merged_pairs/ ./fastqc/ ./linastep_flash2.log ./bowtie/



echo "Finished"
date
