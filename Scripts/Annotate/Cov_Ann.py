import numpy as np
import pandas as pd
from tqdm import tqdm

pos_count2Am = np.load("filtered_cap2_Am.npy", allow_pickle = True)
pos_count1Am = np.load("filtered_cap1_Am.npy", allow_pickle = True)

pos_count2m6Am = np.load("filtered_cap2_m6Am.npy", allow_pickle = True)
pos_count1m6Am = np.load("filtered_cap1_m6Am.npy", allow_pickle = True)

pos_countAm  = np.load("Ambiguous_reads_Am.npy", allow_pickle = True)
pos_countm6Am  = np.load("Ambiguous_reads_m6Am.npy", allow_pickle = True)

tss_seq = pd.read_csv("/Volumes/Extreme SSD/HP/TSNmRNAisoforms.csv", sep = ',')
tss_seq = np.array(tss_seq)

def annotate_to_tss(cap_reads):
    '''
    Annotates reads to transcription start site.

    Input:
        cap_reads (list[list]): list of cap reads to be processed

    Output (array): array of annotated reads
    '''
    annotated = []
    pos_count_dict = {}
    for pos in cap_reads:
        key = (pos[0], pos[1], pos[2], pos[3])
        if key in pos_count_dict:
            pos_count_dict[key][1] += 1  
        else:
            pos_count_dict[key] = [pos, 1] 
    for seq in tqdm(tss_seq):
        key = (seq[0], seq[1], seq[2], seq[3])
        if key in pos_count_dict:
            for i in range(5):
                annotated.append(seq[i])
            annotated.append(pos_count_dict[key][1])
            annotated.append(pos_count_dict[key][0][5])
    annotated = np.array(annotated)
    annotated = annotated.reshape((int(annotated.size/7), 7))
    return annotated

annotated1Am = annotate_to_tss(pos_count1Am)
annotated1m6Am = annotate_to_tss(pos_count1m6Am)

annotated2Am = annotate_to_tss(pos_count2Am)
annotated2m6Am = annotate_to_tss(pos_count2m6Am)

annotatedAm = annotate_to_tss(pos_countAm)
annotatedm6Am = annotate_to_tss(pos_countm6Am)

annotated_both_Am = annotate_to_tss(np.concatenate((pos_count1Am, pos_count2Am)))
annotated_both_m6Am = annotate_to_tss(np.concatenate((pos_count1m6Am, pos_count2m6Am, pos_countm6Am)))

np.save("cap1_ann_Am.npy", annotated1Am)
np.save("cap2_ann_Am.npy", annotated2Am)
np.save("amb_ann_Am.npy", annotatedAm)
np.save("cap12_ann_Am", annotated_both_Am)

np.save("cap1_ann_m6Am.npy", annotated1m6Am)
np.save("cap2_ann_m6Am.npy", annotated2m6Am)
np.save("amb_ann_m6Am.npy", annotatedm6Am)
np.save("cap12_ann_m6Am", annotated_both_m6Am)
