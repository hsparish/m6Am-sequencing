# Description
## Data Analysis workflow:
### Pre-processing

(1) Use cutadapt to remove 3' end DNA adaptor
    
    command used: cutadapt -a AGATCGGAAGAGCGTCGTGT  -A CAGACGTGTGCTCTTCCGATC -m 30 -o R1.cut.fastq -p R2.cut.fastq SRR18074809_1.fastq SRR18074809_2.fastq

(2) Use umi_tools to extract 6 nt UMI tag

    command used: umi_tools extract --bc-pattern=NNNNNN -I R1.cut.fastq --read2-in=R2.cut.fastq -S R1.trim.fastq --read2-out=R2.trim.fastq

(3) Run Eliminating_Duplicates_for_Midway.ipynb script to remove duplicate reads

    (uses barcode pattern and sequence to eliminate duplicate reads, modified version of pyFastDuplicatRemover.py from pyCRAC)

### Labeling reads
(4) Run Finding_caps.ipynb script to label reads as cap1 and cap2

    Script labels reads as cap1 and cap2 using palindrome (allows for 0 or 1 non-template additions)

    Palindrome examples: cap2_1: CCTGAGG, cap2_0: CCTAGG, cap1_1: CGTCG, cap1_0: CGCG

### Aligning reads
(5) Use Hisat-3n to align reads to reference sequence

    Command used: hisat-3n -x hg38 -f -U /scratch/midway2/hparish/Batch-Files/read1_cut.fasta -S cap1_m6am_test.sam --time --base-change A,G --no-spliced-alignment --no-softclip --norc --unique-only --add-chrname --rna-strandness F


### Filtering reads
(6) Run Remove_beginning_mismatch script on cap1 and cap2 bam files

(7) Run separate_ambiguous_cap_assignments script to separate out reads and sites with ambiguous cap1 and cap2 assignments

### Annotating all reads to TSS data
(8) Run cov_ann script to annotate 5' end coverage of reads to TSS-seq data

## Downstream Analysis

(1) Get m6Am stoichiometry info (and percentages in cap1 vs. cap2) with Data_analysis_m6A_m6Am_stoichiometry script

(2) Get global percent methylation at each site relative to the TSN with Global_Percent_Methylation_by_pos script

(3) Annotate m6A sites (and get info about proximity to TSN) with m6A_sites_annotation script

## Important Functions

(1) Function to return base at reference given chromosome and position

(2) Function to separate reads based on gene (assuming TSS data contains this info stored in array)
