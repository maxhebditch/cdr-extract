#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import shutil
import subprocess
import re
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib
import matplotlib.pyplot as plt
from IPython.display import HTML


# In[2]:


results = {}

class antibody:
    def __init__(self, name):
        self.name = name
        self.VH_sequence = VH_sequence
        self.VL_sequence = VL_sequence

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

fasta_file = args.filename

with open(fasta_file,"r") as fasta:
    ab_name = ""
    VH_sequence = ""
    VL_sequence = ""
    VH_found = False
    VL_found = False
    for line in fasta:
        if line.startswith(">"):
            ab_name = line[:-4].replace(">","")
            if "VH" in line:
                VH_found = True
            elif "VL" in line:
                VL_found = True
            else:
                print("not VH or VL")
        else:
            if VH_found == True:
                VH_sequence = line.replace("\n","")
                VH_found = False
            elif VL_found == True:
                VL_sequence = line.replace("\n","")
                VL_found = False
            else:
                print("neither VH nor VL found")
        if ab_name != "" and VH_sequence != "" and VL_sequence != "":
            results[ab_name] = antibody(ab_name)
            ab_name = ""
            VH_sequence = ""
            VL_sequence = ""
        


# In[53]:


def build_CDR_and_non_array(input_sequence,input_index):
    CDR_all_index = []
    for CDR in input_index:
        for i in range(CDR[0],CDR[1]):
            CDR_all_index.append(i)

    CDR_combined_sequence = []
    non_CDR_combined_sequence = []
    for i, aa in enumerate(input_sequence):
        if i in CDR_all_index:
            CDR_combined_sequence.append(aa)
        else:
            non_CDR_combined_sequence.append(aa)
    
    return CDR_combined_sequence, non_CDR_combined_sequence

def find_sub_list(sl,l):
    results=[]
    sll=len(sl)
    for ind in (i for i,e in enumerate(l) if e==sl[0]):
        if l[ind:ind+sll]==sl:
            results.append((ind,ind+sll-1))
    return results

def getindex(cdr,chain):
    indx = find_sub_list(cdr,chain)[0]
    start = indx[0]
    end = indx[1]+1
    return start, end

def CDRH_finder(sequence):
    CDRH_index = []
    aseq = list(sequence)
    #CDRH1
    skip=19
    tseq = aseq[skip:]
    first_C = tseq.index("C")
    four_after_first_C = skip+first_C+4
    start_CDRH1 = four_after_first_C
    rough_CDRH1 = aseq[start_CDRH1:start_CDRH1+13]
    term_W = rough_CDRH1[::-1].index("W")
    CDRH1_a = rough_CDRH1[:-term_W-1]
    pos = find_sub_list(CDRH1_a,aseq)
    if len(pos) > 1:
        print("multiple CDRH1 in sequence")
    else:
        pos = pos[0]
    CDRH1_start = pos[0]
    CDRH1_stop = pos[1]-1
    CDRH1 = "".join(CDRH1_a)
    CDRH1_start, CDRH1_end = getindex(CDRH1_a,aseq)
    if len(CDRH1) > 12 or len(CDRH1) < 10:
        print("CDRH1 not the right size")
        
    #CDRH2
    CDRH2_start = CDRH1_stop + 16
    rough_CDRH2 = aseq[CDRH2_start:CDRH2_start+15]
    end_CDRH2 = rough_CDRH2[-5:]

    for aa in range(0,6):
        close_aa = aseq[CDRH2_start+10+aa]
        far_aa = aseq[CDRH2_start+10+aa+33]
        if far_aa == "C":
            CDRH2_stop = CDRH2_start+10+aa+4
            CDRH3_start = CDRH2_start+10+aa+33+3
    CDRH2_a = aseq[CDRH2_start:CDRH2_stop]
    CDRH2 = "".join(CDRH2_a)
    CDRH2_start, CDRH2_end = getindex(CDRH2_a,aseq)
        
    #CDRH3
    rough_CDRH3 = aseq[CDRH3_start:]
    CDRH3_end_found = False
    multi_pos = find_sub_list(['W','G'],rough_CDRH3)
    if len(multi_pos) > 0:
        for pot_pos in multi_pos:
            if rough_CDRH3[pot_pos[1]+2] == "G":
                CDRH3_end_found = True
                end_CDRH3 = pot_pos[0]
    if CDRH3_end_found == False:
        for ind, aa in enumerate(rough_CDRH3):
            if aa == "W":
                if rough_CDRH3[ind+3] == "G":
                    CDRH3_end_found = True
                    end_CDRH3 = ind
    CDRH3_a = rough_CDRH3[0:end_CDRH3]
    CDRH3 = "".join(CDRH3_a)
    CDRH3_start, CDRH3_end = getindex(CDRH3_a,aseq)
    
    CDRH_index.append([CDRH1_start,CDRH1_end])
    CDRH_index.append([CDRH2_start,CDRH2_end])
    CDRH_index.append([CDRH3_start,CDRH3_end])
    
    
    return CDRH1, CDRH2, CDRH3, CDRH_index
    
def CDRL_finder(sequence):
    CDRL_index = []
    aseq = list(sequence)
    
    #CDRL1
    skip=20
    tseq = aseq[skip:]
    first_C = tseq.index("C")
    CDRL1_start = skip + first_C + 1
    CDRL1_rough = aseq[CDRL1_start:CDRL1_start+20]
    term_W = CDRL1_rough[::-1].index("W")
    CDRL1_a = CDRL1_rough[:-term_W-1]
    CDRL1_start, CDRL1_end = getindex(CDRL1_a,aseq)
    CDRL1 = "".join(CDRL1_a)

    #CDRL2
    CDRL2_start = CDRL1_end + 15
    CDRL2_a = aseq[CDRL2_start:CDRL2_start+7]
    CDRL2_start, CDRL2_end = getindex(CDRL2_a,aseq)
    CDRL2 = "".join(CDRL2_a)
    
    #CDRL3
    CDRL3_start = CDRL2_end + 32
    CDRL3_rough = aseq[CDRL3_start:]
    multi_pos = find_sub_list(['F','G'],CDRL3_rough)
    for pos in multi_pos:
        if CDRL3_rough[pos[1]+2] == "G":
            CDRL3_end_found = True
            CDRL3_end = pos[0]
    if CDRL3_end_found:
        CDRL3_a = aseq[CDRL3_start:CDRL3_start+CDRL3_end]
        CDRL3_start, CDRL3_end = getindex(CDRL3_a,aseq)
        CDRL3 = "".join(CDRL3_a)
        
    else:
        print("CDRL3 end not found")
        
    CDRL_index.append([CDRL1_start,CDRL1_end])
    CDRL_index.append([CDRL2_start,CDRL2_end])
    CDRL_index.append([CDRL3_start,CDRL3_end])
    return CDRL1, CDRL2, CDRL3, CDRL_index
 

for antibody in results:
    VH_sequence = (results[antibody].VH_sequence)
    VL_sequence = (results[antibody].VL_sequence)
    CDRL1, CDRL2, CDRL3, CDRL_index = CDRL_finder(VL_sequence)
    CDRH1, CDRH2, CDRH3, CDRH_index = CDRH_finder(VH_sequence)
    VL_CDR_combined, VL_non_CDR_combined = build_CDR_and_non_array(VL_sequence,CDRL_index)
    VH_CDR_combined, VH_non_CDR_combined = build_CDR_and_non_array(VH_sequence,CDRH_index)
    
    results[antibody].CDRL1 = CDRL1
    results[antibody].CDRL2 = CDRL2
    results[antibody].CDRL3 = CDRL3
    results[antibody].CDRH1 = CDRH1
    results[antibody].CDRH2 = CDRH2
    results[antibody].CDRH3 = CDRH3
    
    results[antibody].VL_CDR_combined = "".join(VL_CDR_combined)
    results[antibody].VL_non_CDR_combined = "".join(VL_non_CDR_combined)
    results[antibody].VH_CDR_combined = "".join(VH_CDR_combined)
    results[antibody].VH_non_CDR_combined = "".join(VH_non_CDR_combined)
    
    results[antibody].all_CDR_combined = results[antibody].VH_CDR_combined+results[antibody].VL_CDR_combined
    results[antibody].all_non_CDR_combined = results[antibody].VH_non_CDR_combined+results[antibody].VL_non_CDR_combined


# In[55]:


with open("all_extracted.fasta","w") as outFile:
    for antibody in results:
        outFile.write(">"+antibody+"_VH\n")
        outFile.write(results[antibody].VH_sequence+"\n")
        outFile.write(">"+antibody+"_VL\n")
        outFile.write(results[antibody].VL_sequence+"\n")
        outFile.write(">"+antibody+"_CDRH1\n")
        outFile.write(results[antibody].CDRH1+"\n")
        outFile.write(">"+antibody+"_CDRH2\n")
        outFile.write(results[antibody].CDRH2+"\n")
        outFile.write(">"+antibody+"_CDRH3\n")
        outFile.write(results[antibody].CDRH3+"\n")
        outFile.write(">"+antibody+"_CDRL1\n")
        outFile.write(results[antibody].CDRL1+"\n")
        outFile.write(">"+antibody+"_CDRL2\n")
        outFile.write(results[antibody].CDRL2+"\n")
        outFile.write(">"+antibody+"_CDRL3\n")
        outFile.write(results[antibody].CDRL3+"\n")
        outFile.write(">"+antibody+"_Heavy_CDRs\n")
        outFile.write(results[antibody].VH_CDR_combined+"\n")
        outFile.write(">"+antibody+"_Light_CDRs\n")
        outFile.write(results[antibody].VL_CDR_combined+"\n")
        outFile.write(">"+antibody+"_Heavy_non_CDRs\n")
        outFile.write(results[antibody].VH_non_CDR_combined+"\n")
        outFile.write(">"+antibody+"_Light_non_CDRs\n")
        outFile.write(results[antibody].VL_non_CDR_combined+"\n")
        outFile.write(">"+antibody+"_Heavy_and_Light_CDRs\n")
        outFile.write(results[antibody].all_CDR_combined+"\n")
        outFile.write(">"+antibody+"_Heavy_and_Light_non_CDRs\n")
        outFile.write(results[antibody].all_non_CDR_combined+"\n")
        
with open("CDRs_separate.fasta","w") as outFile:
    for antibody in results:
        outFile.write(">"+antibody+"_CDRL1\n")
        outFile.write(results[antibody].CDRL1+"\n")
        outFile.write(">"+antibody+"_CDRL2\n")
        outFile.write(results[antibody].CDRL2+"\n")
        outFile.write(">"+antibody+"_CDRL3\n")
        outFile.write(results[antibody].CDRL3+"\n")
        outFile.write(">"+antibody+"_CDRH1\n")
        outFile.write(results[antibody].CDRH1+"\n")
        outFile.write(">"+antibody+"_CDRH2\n")
        outFile.write(results[antibody].CDRH2+"\n")
        outFile.write(">"+antibody+"_CDRH3\n")
        outFile.write(results[antibody].CDRH3+"\n")
        
with open("CDRs_chain_combined.fasta","w") as outFile:
    for antibody in results:
        outFile.write(">"+antibody+"_Heavy_CDRs\n")
        outFile.write(results[antibody].VH_CDR_combined+"\n")
        outFile.write(">"+antibody+"_Light_CDRs\n")
        outFile.write(results[antibody].VL_CDR_combined+"\n")
        
with open("non_CDRs_chain_combined.fasta","w") as outFile:
    for antibody in results:
        outFile.write(">"+antibody+"_Heavy_CDRs\n")
        outFile.write(results[antibody].VH_non_CDR_combined+"\n")
        outFile.write(">"+antibody+"_Light_CDRs\n")
        outFile.write(results[antibody].VL_non_CDR_combined+"\n")
        
with open("CDRs_combined.fasta","w") as outFile:
    for antibody in results:
        outFile.write(">"+antibody+"_Heavy_and_Light_CDRs\n")
        outFile.write(results[antibody].all_CDR_combined+"\n")
        
with open("non_CDRs_combined.fasta","w") as outFile:
    for antibody in results:
        outFile.write(">"+antibody+"_Heavy_and_Light_non_CDRs\n")
        outFile.write(results[antibody].all_non_CDR_combined+"\n")

