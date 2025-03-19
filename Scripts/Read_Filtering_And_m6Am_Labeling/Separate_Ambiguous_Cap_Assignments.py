import numpy as np
from tqdm import tqdm

cap1Am = np.load("cap1_Am_nm.npy", allow_pickle = True)
cap1m6Am = np.load("cap1_m6Am_nm.npy", allow_pickle = True)

cap2Am = np.load("cap2_Am_nm.npy", allow_pickle = True)
cap2m6Am = np.load("cap2_m6Am_nm.npy", allow_pickle = True)

pos_count2m6Am_filt = []
pos_count1Am_filt = []
pos_count1m6Am_filt = []

cap1 = np.concatenate((cap1Am, cap1m6Am))
cap2 = np.concatenate((cap2Am, cap2m6Am))

def pos_dict(sequences):
    '''
    Constructs a dictionary of the positions of the reads in a list of sequences.

    Input:
        sequences (list[list]): list of sequence information for reads
    Output (dict): dictionary of read positions with sequence as the key and a tuple 
        with chromosome name, position, and strand as the value
    '''
    pos_count_dict = {}
    for pos in sequences:
        key = ((pos[4])[:40])
        pos_count_dict[key] = pos[0:5]
    return pos_count_dict

def find_ambiguous(sequences, pos_count_dict):
    '''
    Finds the reads that were labeled as both cap1 and cap2.

    Input:
        sequences (list[list]): list of sequence information for reads
        pos_count_dict (dict): dictionary of read positions
    Output (list, dict): list of ambiguous reads and dictionary of ambiguous reads
    '''
    for seq in tqdm(sequences):
        pos_count_dups = []
        key = ((seq[4])[:40])
        if key in pos_count_dict:
            pos_count_dups.append(seq[0:5])
    pos_count_dups_dict = {}
    for pos in pos_count_dups:
        key = tuple(pos[0:4])
        pos_count_dups_dict[key] = pos[0:5]
    return pos_count_dups, pos_count_dups_dict

def ambiguous_positions(sequences, pos_count_dups_dict, pos_count_dups):
    '''
    Constructs dictionary with the positions of the ambiguous reads as key.

    Input:
        sequences (list[list]): list of sequence information for reads
        pos_count_dups_dict (dict): dictionary of ambiguous reads
        pos_count_dups (list): list of ambiguous reads
    Output (dict): dictionary of ambiguous reads and any other reads matching those positions
    '''
    pos_count_dups_remove = {}
    for seq in sequences:
        key = tuple(seq[0:4])
        if key in pos_count_dups_dict:
            pos_count_dups.append(seq[0:5])
    for pos in pos_count_dups:
        key = tuple(pos[0:5])
        pos_count_dups_remove[key] = pos[0:5]
    return pos_count_dups_remove

def filter_reads(sequences, pos_count_dups_remove):
    '''
    Removes ambiguous reads and reads matching those positions.

    Input:
        sequences (list[list]): list of sequence information for reads
        pos_count_dups (list): list of ambiguous reads to remove
    Output (list, list): list of filtered reads and list of ambiguous reads
    '''
    filtered = []
    ambiguous = []
    for seq in tqdm(sequences):
        key = tuple(seq[0:5])
        if key not in pos_count_dups_remove:
            filtered.append(seq)
        if key in pos_count_dups_remove:
            ambiguous.append(seq)
    return filtered, ambiguous

pos_count2_dict = pos_dict(cap2) 
pos_count1_dict = pos_dict(cap1) 
pos_count2_dups, pos_count2_dups_dict = find_ambiguous(cap2, pos_count1_dict)
pos_count1_dups, pos_count1_dups_dict = find_ambiguous(cap1, pos_count2_dict)
pos_count2_dups_remove = ambiguous_positions(cap2, pos_count1_dups_dict, pos_count2_dups)
pos_count1_dups_remove = ambiguous_positions(cap1, pos_count2_dups_dict, pos_count1_dups)
pos_count2Am_filt, ambAm = filter_reads(cap2Am, pos_count2_dups_remove)
pos_count2m6Am_filt, ambm6Am = filter_reads(cap2m6Am, pos_count2_dups_remove)
pos_count1Am_filt, _ = filter_reads(cap1Am, pos_count1_dups_remove)
pos_count1m6Am_filt, _ = filter_reads(cap1m6Am, pos_count1_dups_remove)

pos_count2Am_filt = np.array(pos_count2Am_filt)
pos_count1Am_filt = np.array(pos_count1Am_filt)
pos_count2m6Am_filt = np.array(pos_count2m6Am_filt)
pos_count1m6Am_filt = np.array(pos_count1m6Am_filt)
ambAm = np.array(ambAm)
ambm6Am = np.array(ambm6Am)

np.save("filtered_cap1_Am.npy", pos_count1Am_filt)
np.save("filtered_cap2_Am.npy", pos_count2Am_filt)
np.save("filtered_cap1_m6Am.npy", pos_count1m6Am_filt)
np.save("filtered_cap2_m6Am.npy", pos_count2m6Am_filt)
np.save("ambiguous_reads_Am.npy", ambAm)
np.save("ambiguous_reads_m6Am.npy", ambm6Am)
