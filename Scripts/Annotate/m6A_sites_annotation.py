import pysam
import numpy as np
import pandas as pd
from tqdm import tqdm
from pyfaidx import Fasta
from ipynb.fs.full. Important_Functions import *

Am2 = np.load("cap2_ann_Am.npy", allow_pickle = True)
Am1 = np.load("cap1_ann_Am.npy", allow_pickle = True)

m6Am2 = np.load("cap2_ann_m6Am.npy", allow_pickle = True)
m6Am1 = np.load("cap1_ann_m6Am.npy", allow_pickle = True)

Am  = np.load("amb_ann_Am.npy", allow_pickle = True)
m6Am  = np.load("amb_ann_m6Am.npy", allow_pickle = True)

file_path = "/Volumes/Extreme SSD/HP/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
fasta = Fasta(file_path)

def A_m6A_array(tsn_mod_array):
    '''
    Locate sites of m6A modifications in the reads.
    
    Input:
        tsn_mod_array (list): list of annotated reads
    Output ((list, list)): tuple with list of m6A sites and list of A sites
    '''
    m6A_array = []
    A_array = []
    for i in range(len(tsn_mod_array)):
        for j in range(len(tsn_mod_array[i][6])):
            if j > 0 and tsn_mod_array[i][6][j] == 'A':
                m6A_array.append(tsn_mod_array[i][0])
                m6A_array.append(int(tsn_mod_array[i][1])+j+1)
                m6A_array.append(int(tsn_mod_array[i][2])+j+1)
                m6A_array.append(tsn_mod_array[i][3])
                m6A_array.append(tsn_mod_array[i][4])
                m6A_array.append(j+1)                  
            if j > 0 and tsn_mod_array[i][6][j] == 'G' and chr_loc(str(tsn_mod_array[i][0]), int(tsn_mod_array[i][1]), fasta) == 'A':
                A_array.append(tsn_mod_array[i][0])
                A_array.append(int(tsn_mod_array[i][1])+j+1)
                A_array.append(int(tsn_mod_array[i][2])+j+1)
                A_array.append(tsn_mod_array[i][3])
                A_array.append(tsn_mod_array[i][4])
                A_array.append(j+1)  
    m6A_array = np.array(m6A_array)
    m6A_array = m6A_array.reshape((int(m6A_array.size/6), 6))
    
    A_array = np.array(A_array)
    A_array = A_array.reshape((int(A_array.size/6), 6))
    return m6A_array, A_array

m6A_Am2, A_Am2 = A_m6A_array(Am2)

m6A_m6Am2, A_m6Am2 = A_m6A_array(m6Am2) 

m6A_Am1, A_Am1 = A_m6A_array(Am1) 

m6A_m6Am1, A_m6Am1 = A_m6A_array(m6Am1) 

def annotate(mod_array):
    '''
    Find number of reads for each m6A modification vs. A modification.

    Input:
        mod_array (list): list of reads with A/m6A sites labeled
    Output (array): sites annotated with following info:
        chromosome
    '''
    dict_sites = {}
    for elt in mod_array:
        key = (elt[0], elt[1], elt[2], elt[3], elt[5])
        if key in dict_sites:
            dict_sites[key][1] += 1  
        else:
            dict_sites[key] = [elt, 1] 

    mod_sites = []
    mod_ann = [(k, v) for k, v in dict_sites.items()]
    for i in range(len(mod_ann)):
        for j in range(6):
            mod_sites.append(mod_ann[i][1][0][j])
        mod_sites.append(mod_ann[i][1][1])

    mod_sites = np.array(mod_sites)
    mod_sites = mod_sites.reshape((int(mod_sites.size/7), 7))
    return mod_sites

m6A_Am1_sites = annotate(m6A_Am1)
A_Am1_sites = annotate(A_Am1)
m6A_m6Am1_sites = annotate(m6A_m6Am1)
A_m6Am1_sites = annotate(A_m6Am1)
m6A_Am2_sites = annotate(m6A_Am2)
A_Am2_sites = annotate(A_Am2)
m6A_m6Am2_sites = annotate(m6A_m6Am2)
A_m6Am2_sites = annotate(A_m6Am2)

np.save("A_m6Am1_ann.npy", A_m6Am1_sites)
np.save("m6A_m6Am1_ann.npy", m6A_m6Am1_sites)
np.save("A_Am1_ann.npy", A_Am1_sites)
np.save("m6A_Am1_ann.npy", m6A_Am1_sites)

np.save("A_m6Am2_ann.npy", A_m6Am2_sites)
np.save("m6A_m6Am2_ann.npy", m6A_m6Am2_sites)
np.save("A_Am2_ann.npy", A_Am2_sites)
np.save("m6A_Am2_ann.npy", m6A_Am2_sites)
