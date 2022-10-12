### 1.1 Introduction to the command line for Bioinformatics and Quality control with the FASTQC tool ###
# try basic Linux/OSX commands, such as
ls
cd data2022

# observe the first ten lines -what do you notice?
head Platz10_R1.head.fastq

# count the total number of lines with the command below -how many READS are in this FASTQ?
# hint: each read in FASTQ format consists of four lines
wc -l Platz10_R1.head.fastq

# find specific words, e.g "AATATT"
grep "AATATT" Platz10_R1.head.fastq

# combine commands. e.g.
grep "AATATT" Platz10_R1.head.fastq | wc -l

# quality check with fastqc
fastqc --quiet Platz1_R1.head.fastq
# click to HTML output and select "View in Web Browser"
# a new tab opens with the fastqc results
# repeat the analyses with the second sequencing pair, e.g. Platz1_R2.head.fastq

### 1.2 Amplicon assembly with pandaseq ###

# pandaseq options
pandaseq -h

# amplicon assembly with pandaseq
pandaseq -f Platz10_R1.head.fastq -r Platz10_R2.head.fastq -w Platz10.fa -g log.txt

# observe the first ten lines of the FASTA output
head Platz10.fa

# how many sequences are in this FASTA?
# hint: count the total number of HEADER line symbols (">") with any of the commands below 
grep ">" Platz10.fa | wc -l
grep -c ">" Platz10.fa

# Q: what is the rate of FASTA sequences vs FASTQ reads?
# and what does this tell us about our sequencing and assembly quality and efficiency?

### 1.3 OTU Assignment with kraken2 against the 16S Greengenes database ###

# Current path: ~/data2022
# Note: Run this step together with your neighbor/colleague to avoid saturation of the storage server

# Run your samples against Greengenes with your own data
# the example here is the with the Platz10 run
kraken2 --db 16S_Greengenes_k2db --use-names --output output.txt --report report.txt Platz10.fa

# kraken can also run directly with FASTQ files
# replace the FASTQ filenames with your own FASTQ names
kraken2 --db 16S_Greengenes_k2db --use-names --output output.txt --report report.txt --paired Platz10_R1.head.fastq Platz10_R2.head.fastq

# inspect your individual results (file: report.txt) within the RStudio window
# Q: What are the most predominant genera in your personal Illumina runs?
