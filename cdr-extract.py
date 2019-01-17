#!/usr/bin/env python

import os
import re
import argparse

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
    
    return "".join(CDR_combined_sequence), "".join(non_CDR_combined_sequence)

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
 
all_VL = {}
all_VH = {}
all_chains = {}

class format_result:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

fasta_file = args.filename
fasta_file_name = os.path.splitext(fasta_file)[0]

sequence_found = False
with open(fasta_file,"r") as fasta:
    for line in fasta:
        #Look for name line
        if line.startswith(">"):
            if sequence_found == True:
                sequence = "".join(frag for frag in chain_array)
                all_chains[ab_name] = format_result(ab_name,sequence)
                sequence_found = False
                chain_array = []
            else:
                chain_array = []
            ab_name = line.replace(">","")
            ab_name = ab_name.strip()
        #If not a name line
        else:
            sequence_found = True
            chain_sequence = line.strip()
            chain_array.append(chain_sequence)

if len(chain_array) > 0:
    sequence = "".join(frag for frag in chain_array)
    all_chains[ab_name] = format_result(ab_name,sequence)
    sequence_found = False
    chain_array = []

for result in all_chains:
    ab_name = all_chains[result].name
    print("Processing "+str(ab_name))
    chain_str = (all_chains[result].sequence)

    heavy_pattern = re.compile("VS[SA]")
    light_pattern = re.compile("(E[LIV]KR)|(E[LIV]KKR)|EIIKR|(TVL[GSA])|(ERKR)")  

    def extract_V(match_array):
        if len(match_array) > 1:
            V_string = "".join(seq for seq in match_array)
        else:
            V_string = match_array[0]
        
        return V_string

    if re.search(light_pattern, chain_str) is not None:
        match = re.search(light_pattern, chain_str).group(0)
        light_regex = r"(.*?)" + re.escape(match)
        VL_array = re.findall(light_regex, chain_str)
        VL_sequence = extract_V(VL_array)
        all_VL[ab_name] = format_result(ab_name,VL_sequence)
        print("Found light chain")
    elif re.search(heavy_pattern, chain_str) is not None:
        VH_array = re.findall("(.*?)VS[SA]", chain_str)
        VH_sequence = extract_V(VH_array)
        all_VH[ab_name] = format_result(ab_name,VH_sequence)
        print("Found heavy chain")
    else:
        print("Suspected non-ab chain")

if len(all_VL) > 0:
    with open(fasta_file_name+"_VL_CDR"+".fasta","w") as outFile:
        print("")
        for VL in all_VL:
            print("Writing light chain "+str(all_VL[VL].name))
            VL_sequence = all_VL[VL].sequence
            CDRL1, CDRL2, CDRL3, CDRL_index = CDRL_finder(VL_sequence)
            VL_CDR_combined, VL_non_CDR_combined = build_CDR_and_non_array(VL_sequence,CDRL_index)
            
            outFile.write(">"+all_VL[VL].name+"_VL\n")
            outFile.write(all_VL[VL].sequence+"\n")
            outFile.write(">"+all_VL[VL].name+"_CDRL1\n")
            outFile.write(CDRL1+"\n")
            outFile.write(">"+all_VL[VL].name+"_CDRL2\n")
            outFile.write(CDRL2+"\n")
            outFile.write(">"+all_VL[VL].name+"_CDRL3\n")
            outFile.write(CDRL3+"\n")
            outFile.write(">"+all_VL[VL].name+"_ALL_L_CDR\n")
            outFile.write(VL_CDR_combined+"\n")
            outFile.write(">"+all_VL[VL].name+"_non_L_CDR\n")
            outFile.write(VL_non_CDR_combined+"\n")

if len(all_VH) > 0:
    with open(fasta_file_name+"_VH_CDR"+".fasta","w") as outFile:
        print("")
        for VH in all_VH:
            print("Writing heavy chain "+str(all_VH[VH].name))
            VH_sequence = all_VH[VH].sequence
            CDRH1, CDRH2, CDRH3, CDRH_index = CDRH_finder(VH_sequence)
            VH_CDR_combined, VH_non_CDR_combined = build_CDR_and_non_array(VH_sequence,CDRH_index)

            outFile.write(">"+all_VH[VH].name+"_VH\n")
            outFile.write(all_VH[VH].sequence+"\n")
            outFile.write(">"+all_VH[VH].name+"_CDRH1\n")
            outFile.write(CDRH1+"\n")
            outFile.write(">"+all_VH[VH].name+"_CDRH2\n")
            outFile.write(CDRH2+"\n")
            outFile.write(">"+all_VH[VH].name+"_CDRH3\n")
            outFile.write(CDRH3+"\n")
            outFile.write(">"+all_VH[VH].name+"_combined_H_CDR\n")
            outFile.write(VH_CDR_combined+"\n")
            outFile.write(">"+all_VH[VH].name+"_combined_non_H_CDR\n")
            outFile.write(VH_non_CDR_combined+"\n")
