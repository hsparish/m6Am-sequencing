##Data Analysis workflow:

(1) Use cutadapt to remove 3' end DNA adaptor

(2) Use umi_tools to extract 6 nt UMI tag

(3) Run Eliminating_Duplicates_for_Midway.ipynb script to remove duplicate reads

(4) Run Finding_caps.ipynb script to label reads as cap1 and cap2

(5) Use Hisat-3n to align reads to reference sequence

(6) 
