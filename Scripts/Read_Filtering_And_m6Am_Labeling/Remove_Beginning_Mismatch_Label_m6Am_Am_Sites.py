import pysam
import numpy as np
from tqdm import tqdm
from pyfaidx import Fasta
from ipynb.fs.full. Important_Functions import *

samfile1 = pysam.AlignmentFile("/Volumes/Extreme SSD/HP/cap1_m6am_test.bam", "rb")
samfile2 = pysam.AlignmentFile("/Volumes/Extreme SSD/HP/cap2_m6am_test.bam", "rb")
ref_seq = "/Volumes/Extreme SSD/HP/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
reference = Fasta(ref_seq)

def fetch_reads(samfile):
    '''
    Retrieve reads from provided samfile.

    Input:
        samfile (str): file path to samfile
    Output (List[List[reads]]): read information from samfile by chromosome
    '''
    reads = []
    for i in range(24):
        reads.append([])
    for i in range(1,25):
        if i == 23:
            for read in samfile.fetch('chrX'):
                reads[i-1].append(read) 
        elif i == 24:
            for read in samfile.fetch('chrY'):
                reads[i-1].append(read) 
        else:
            for read in samfile.fetch(f'chr{i}'):
                reads[i-1].append(read)
    return reads

def add_read_info(read, pos_count_array_m6Am, pos_count_array_Am, pos_count_array, i, direction):                    
    '''
    Modifies provided array by appending the desired read information.

    Input:
        read: samfile read
        pos_count_array_m6Am (list): list to which m6Am reads are appended        
        pos_count_array_Am (list): list to which Am reads are appended
        pos_count_array (list): list to which all other reads are appended
        i (int): chromosome number
        direction (str): forward or reverse strand
    '''
    if direction == 'fwd':
        n = 'A'
        ref = 'G'
    else:
        n = 'T'
        ref = 'C'
    if read.query_sequence[0] == n:
        format_reads(read, pos_count_array_m6Am, i, 'fwd')
    elif read.query_sequence[0] == ref and chr_loc(f'chr{i+1}', read.reference_end-39, reference) == n:
        format_reads(read, pos_count_array_Am, i, 'fwd')
    else:
        format_reads(read, pos_count_array, i, 'fwd')           

def format_reads(read, formatted, i, direction):
    '''
    Formats arrays containing the desired read information.

    Input:
        read: samfile read
        formatted (list): list to which read information is being appended
        i (int): chromosome number
        direction (str): forward or reverse strand
    '''
    if i == 22:
        formatted.append('chrX')
    elif i == 23:
        formatted.append('chrY')
    else:
        formatted.append(f'chr{i+1}')
    if direction == 'fwd':
        formatted.append(read.reference_end-40)
        formatted.append(read.reference_end-39)
        formatted.append("+")
    elif direction == 'rev':
        formatted.append(read.reference_end-1)
        formatted.append(read.reference_end)
        formatted.append("-")
    formatted.append(str(read.query_name).split("_")[-1] + read.query_sequence)   
    formatted.append(str(read.query_sequence))   

def check_sequence(read, ref_chr, ref_pos, num_bases, is_forward):
    '''
    Check if the read's sequence matches the reference sequence at the given position.
    
    Input:
        read: samfile read
        ref_chr (str): name of chromosome
        ref_pos (int): position along reference to check
        num_bases (int): 3 for cap2, 2 for cap1
        is_forward (boolean): True if it is a forward sequence
    Ouptut (boolean): True if there is no mismatch at the beginning of read
    '''
    for k in range(num_bases):
        read_base = read.query_sequence[k] if is_forward else read.query_sequence[39 - k]
        ref_base = chr_loc(ref_chr, ref_pos + k if is_forward else ref_pos - k, reference)
        if not (read_base == ref_base or (read_base, ref_base) in [('G', 'A'), ('C', 'T')]):
            return False
    return True

def process_reads(reads, pos_count_m6Am, pos_count_Am, pos_count, i, num_bases):
    '''
    Process reads for a given chromosome index and number of bases to compare.

    Input:
        reads (List): list of samfile reads
        pos_count_m6Am (List): list of m6Am reads to be modified
        pos_count_Am (List): list of Am reads to be modified
        pos_count (List): list of all other reads to be modified
        i (int): chromosome number
        num_bases (int): 3 for cap2, 2 for cap1
    '''  
    ref_chr = f'chr{i+1}'
    for read in reads:
        ref_pos = read.reference_end - (39 if read.is_forward else 0)
        if check_sequence(read, ref_chr, ref_pos, num_bases, read.is_forward):
            strand = "fwd" if read.is_forward else "rev"
            add_read_info(read, pos_count_m6Am, pos_count_Am, pos_count, i, strand)

reads1 = fetch_reads(samfile1)
reads2 = fetch_reads(samfile2)

pos_count1, pos_count1_m6Am, pos_count1_Am = [], [], []
pos_count2, pos_count2_m6Am, pos_count2_Am = [], [], []

for i in tqdm(range(24)):
    process_reads(reads1[i], pos_count1_m6Am, pos_count1_Am, pos_count1, i, num_bases=2)
    process_reads(reads2[i], pos_count2_m6Am, pos_count2_Am, pos_count2, i, num_bases=3)

cap1Am, cap1m6Am = np.array(pos_count1_Am, dtype = 'object')
cap1m6Am = np.array(pos_count1_m6Am, dtype = 'object')

cap2Am = np.array(pos_count2_Am, dtype = 'object')
cap2m6Am = np.array(pos_count2_m6Am, dtype = 'object')
      
cap1Am = cap1Am.reshape((int(cap1Am.size/6), 6))
cap2Am = cap2Am.reshape((int(cap2Am.size/6), 6))

cap1m6Am = cap1m6Am.reshape((int(cap1m6Am.size/6), 6))
cap2m6Am = cap2m6Am.reshape((int(cap2m6Am.size/6), 6))

np.save("cap1_Am_nm.npy", cap1Am)
np.save("cap2_Am_nm.npy", cap2Am)
np.save("cap1_m6Am_nm.npy", cap1m6Am)
np.save("cap2_m6Am_nm.npy", cap2m6Am)