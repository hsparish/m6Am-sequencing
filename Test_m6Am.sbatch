#!/bin/bash

#SBATCH --job-name=single-node-cpu-example
#SBATCH --account=pi-weixintang
#SBATCH --partition=bigmem2  # accessible partitions listed by the sinfo command
#SBATCH --ntasks-per-node=1  # number of tasks per node
#SBATCH --cpus-per-task=16    # number of CPU cores per task

# Load the require module(s)
module load samtools

# match the no. of threads with the no. of CPU cores
/home/hparish/hisat-3n/hisat-3n -x hg38 -f -U /scratch/midway2/hparish/Batch-Files/read1_cut.fasta -S cap1_m6am_test.sam --time --base-change A,G --no-spliced-alignment --no-softclip --norc --no-unal --add-chrname --rna-strandness F

samtools sort cap1_m6am_test.sam -o cap1_m6am_test.bam

samtools index cap1_m6am_test.bam



/home/hparish/hisat-3n/hisat-3n -x hg38 -f -U /scratch/midway2/hparish/Batch-Files/read_2.fasta -S cap2_m6am_test.sam --time --base-change A,G --no-spliced-alignment --no-softclip --norc --no-unal --add-chrname --rna-strandness F

samtools sort cap2_m6am_test.sam -o cap2_m6am_test.bam

samtools index cap2_m6am_test.bam

