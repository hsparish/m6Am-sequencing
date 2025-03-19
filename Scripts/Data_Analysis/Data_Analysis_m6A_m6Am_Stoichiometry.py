import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

cap1_Am = np.load("cap1_ann_Am.npy", allow_pickle = True)
cap1_m6Am = np.load("cap1_ann_m6Am.npy", allow_pickle = True)
m6Am = np.load("cap12_ann_m6Am.npy", allow_pickle = True)
cap2_Am = np.load("cap2_ann_Am.npy", allow_pickle = True)
cap2_m6Am = np.load("cap2_ann_m6Am.npy", allow_pickle = True)
Am = np.load("cap12_ann_Am.npy", allow_pickle = True)

def find_percents(m6Am_array, Am_array):
    '''
    Calculate percent Am vs m6Am at each TSS.

    Input:
        m6Am_array (list): list of m6Am reads
        Am_array (list): list of Am reads
    Output:
    '''
    overlap = []
    only_Am = []
    only_m6Am = []
    Am_dict = {tuple(Am_array[i][:4]): Am_array[i][5] for i in range(len(Am_array))}
    m6Am_dict = {tuple(m6Am_array[i][:4]): m6Am_array[i][5] for i in range(len(m6Am_array))}
    for i in tqdm(range(len(m6Am_array))):
        key = tuple(m6Am_array[i][:4])
        if key in Am_dict:
            for j in range(6):
                overlap.append((m6Am_array[i])[j])
            overlap.append(Am_dict[key])
            overlap.append(100*int(m6Am_array[i][5])/(int(m6Am_array[i][5])+int(Am_dict[key])))
        else:
            for j in range(6):
                only_m6Am.append((m6Am_array[i])[j])
            only_m6Am.append(0)
            only_m6Am.append(100)

    for i in tqdm(range(len(Am_array))):
        key = tuple(Am_array[i][:4])
        if key in m6Am_dict:
            pass
        else: 
            for j in range(5):
                only_Am.append((Am_array[i])[j])
            only_Am.append(0)
            only_Am.append(Am_array[i][5])
            only_Am.append(0)
            
    overlap = np.array(overlap)
    overlap = overlap.reshape((int(overlap.size/8),8))
    
    only_Am = np.array(only_Am)
    only_Am = only_Am.reshape((int(only_Am.size/8),8))
    
    only_m6Am = np.array(only_m6Am)
    only_m6Am = only_m6Am.reshape((int(only_m6Am.size/8),8))

    return overlap, only_Am, only_m6Am

overlap1, only_Am1, only_m6Am1 = find_percents(cap1_m6Am, cap1_Am)
overlap2, only_Am2, only_m6Am2 = find_percents(cap2_m6Am, cap2_Am)
overlap12, only_Am12, only_m6Am12 = find_percents(m6Am, Am)



