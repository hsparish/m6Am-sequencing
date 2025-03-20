import os
from tqdm import tqdm
import numpy as np
from Bio import SeqIO

def read_fasta(file_path): 
    '''
    Takes a fasta file as input and outputs an array of sequences in the fasta file.
    Input:
        Fasta file path (str)
    Output (list, list):
        Tuple with in which the first item is a list of sequences in the fasta file 
        and the second is a list of sequence IDs
    '''
    sequences = []
    names = []
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file_path, "fasta"):
            names.append(record.id)
            sequences.append(record.seq)
    return sequences, names

file_path = "/Volumes/Extreme SSD/HP/Scripts/dedup+_1.fasta"
sequences, names = read_fasta(file_path)

def reverse_complement(seq):
    '''
    Generates the reverse complement of a given input sequence.

    Input:
        seq (str): DNA sequence
    Output (str): reverse complement of the sequence
    '''
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[n] for n in seq)

def generate_palindromes(n, maximum):
    '''
    Generate the palindrome sequences for cap1 and cap2 structures.

    Input:
        n (int): 3 for cap2, 2 for cap1
        maximum (int): 3 for cap2, 1 for cap1
    Output:
        palindromes (set): list of palindrome sequences for cap1 or cap2
    '''
    if n == 1:
        return ['AT', 'GC', 'TA', 'CG']
    else:
        palindromes = set()
        for x in ['A', 'G', 'T', 'C']:
            for y in generate_pal(n-1, maximum):
                candidate = x + y + reverse_complement(x)
                if len(set(candidate[:maximum])) > 1:  
                    palindromes.add(candidate)  
        return palindromes
    
palindromes_cap2 = generate_palindromes(3)
palindromes_cap1 = generate_palindromes(2)

def match_palindrome(seqs, palindromes, names, i):
    '''
    Finds sequences with a beginning match to cap1 or cap palindromes with 0 or 1 
    nontemplate additions

    Input:
        seqs (List[Sequences]): list of DNA reads
        palindromes (List[str]): list of DNA palindromes
        names (List[str]): list of IDs for reads
        i (int): 2 for cap1, 3 for cap2
    Output ((List[Sequence], List[str], List[Sequence], List[str])): 
        List of sequences with 1 nontemplate addition, list of IDs for those reads, list of sequences with 0 nontemplate additions, lsi
    '''
    sequence1 = []
    sequence0 = []
    name1 = []
    name0 = []
    for j, seq in enumerate(seqs):    # i = 3 cap2 
        match1 = seq[:i] + seq[i+1:2*i+1]
        match0 = seq[:2*i]
        if match1 in palindromes:
            sequence1.append(seq)
            name1.append(names[j])
        elif match0 in palindromes:
            sequence0.append(seq)
            name0.append(names[j])
    return sequence1, name1, sequence0, name0 

cap2_1, cap2_1_name, cap2_0, cap2_0_name = match_palindrome(sequences, palindromes_cap2, names, 3)
cap1_1, cap1_1_name, cap1_0, cap1_0_name = match_palindrome(sequences, palindromes_cap1, names, 2)

def write_fasta(seq0, seq1, name0, name1, file_name, i):
    '''
    Writes labeled sequences to a new fasta file with specified file name.

    Input:
        seq0 (List[Sequence]): list of sequences with no nontemplate additions
        seq1 (List[Sequence]): list of sequences with one nontemplate addition
        name0 (List[str]): list of IDs for sequences with no nontemplate additions
        name1 (List[str]): list of IDs for sequences with one nontemplate additions
        i (int): 2 if cap1, 3 if cap2
    '''
    with open(file_name, 'w') as file:
        for j, seq in tqdm(enumerate(seq0)):
            if len(seq)>=40 + i: #42
                file.write(f">Sequence{j+1}_{name0[j]}_{seq[:i]}\n")
                file.write(str(seq[i:i+40] + "\n"))
        for j, seq in tqdm(enumerate(seq1)):
            if len(seq) >= (40 + i + 1):
                file.write(f">Sequence{j+1}_{name1[j]}_{seq[:i+1]}\n")
                file.write(str(seq[i+1:40+i+1] + "\n"))

write_fasta(cap1_0, cap1_1, cap1_0_name, cap1_1_name,\
            "/Volumes/Extreme SSD/HP/cap1_m6am_test_cut+.fasta", 2)
write_fasta(cap2_0, cap2_1, cap2_0_name, cap2_1_name,\
            "/Volumes/Extreme SSD/HP/cap2_m6am_test_cut+.fasta", 3)

